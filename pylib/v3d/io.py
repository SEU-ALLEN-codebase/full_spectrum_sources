import os
import struct
import numpy as np
import sys


def load_v3draw(path):
    """
    by Zuohan Zhao
    from basic_c_fun/stackutils.cpp
    2022/5/8
    """
    assert os.path.exists(path)
    formatkey = "raw_image_stack_by_hpeng"
    with open(path, "rb") as f:
        filesize = os.path.getsize(path)
        header_len = len(formatkey) + 2 + 4 * 4 + 1
        assert (filesize >= header_len)
        format = f.read(len(formatkey)).decode('utf-8')
        assert (format == formatkey)
        endianCodeData = f.read(1).decode('utf-8')
        if endianCodeData == 'B':
            endian = '>'
        elif endianCodeData == 'L':
            endian = '<'
        else:
            raise Exception('endian be big/little')
        datatype = struct.unpack(endian + 'h', f.read(2))[0]
        if datatype == 1:
            dt = 'u1'
        elif datatype == 2:
            dt = 'u2'
        elif datatype == 4:
            dt = 'f4'
        else:
            raise Exception('datatype be 1/2/4')
        sz = struct.unpack(endian + 'iiii', f.read(4 * 4))
        tot = sz[0] * sz[1] * sz[2] * sz[3]
        tot_buffer_len = tot * datatype
        tot_len = tot_buffer_len + header_len
        if tot_len != filesize:
            f.seek(-4 * 2, 1)
            tot = sz[0] * sz[1] * sz[2] * sz[3]
            assert (tot_len == filesize)
        img = np.frombuffer(f.read(tot_buffer_len), endian + dt)
        return img.reshape(sz[-1:-5:-1])


def save_v3draw(img: np.ndarray, path):
    """
    by Zuohan Zhao
    from basic_c_fun/stackutils.cpp
    2022/5/8
    """
    with open(path, 'wb') as f:
        formatkey = "raw_image_stack_by_hpeng"
        f.write(formatkey.encode())
        if img.dtype.byteorder == '>':
            endian = 'B'
        elif img.dtype.byteorder == '<':
            endian = 'L'
        elif img.dtype.byteorder == '|':
            endian = 'B'
        else:
            if sys.byteorder == 'little':
                endian = 'L'
            else:
                endian = 'B'
        f.write(endian.encode())
        if img.dtype == np.uint8:
            datatype = 1
        elif img.dtype == np.uint16:
            datatype = 2
        else:
            datatype = 4
        endian = '>' if endian == 'B' else '<'
        f.write(struct.pack(endian + 'h', datatype))
        sz = list(img.shape)
        sz.extend([1] * (4 - len(sz)))
        sz.reverse()
        f.write(struct.pack(endian + 'iiii', *sz))
        f.write(img.tobytes())


class PBD:
    """
    by Zuohan Zhao
    from neuron_annotator/utility/ImageLoaderBasic.cpp
    2022/6/23
    """

    def __init__(self):
        self.total_read_bytes = 0
        self.max_decompression_size = 0
        self.channel_len = 0
        self.compression_buffer = bytearray()
        self.decompression_buffer = bytearray()
        self.compression_pos = 0
        self.decompression_pos = 0
        self.decompression_prior = 0
        self.pbd3_src_min = 0
        self.pbd3_src_max = 0
        self.pbd3_cur_min = 0
        self.pbd3_cur_max = 0
        self.pbd3_cur_chan = 0
        self.load_datatype = 0
        self.pbd_sz = [0, 0, 0, 0]
        self.endian = ""

    def decompress_pbd8(self, compression_len: int):
        cp = 0
        dp = 0
        pva = 0
        pvb = 0
        mask = 3
        while cp < compression_len:
            value = self.compression_buffer[self.compression_pos + cp]
            if value < 33:
                count = value + 1
                self.decompression_buffer[self.decompression_pos + dp:self.decompression_pos + dp + count] = \
                    self.compression_buffer[self.compression_pos + cp + 1:self.compression_pos + cp + 1 + count]
                cp += count + 1
                dp += count
                self.decompression_prior = self.decompression_buffer[self.decompression_pos + dp - 1]
            elif value < 128:
                left_to_fill = value - 32
                while left_to_fill > 0:
                    fill_num = min(left_to_fill, 4)
                    cp += 1
                    src_char = self.compression_buffer[self.compression_pos + cp]
                    to_fill = self.decompression_pos + dp
                    p0 = src_char & mask
                    src_char >>= 2
                    p1 = src_char & mask
                    src_char >>= 2
                    p2 = src_char & mask
                    src_char >>= 2
                    p3 = src_char & mask
                    pva = self.decompression_prior + (-1 if p0 == 3 else p0)
                    self.decompression_buffer[to_fill] = pva
                    if fill_num > 1:
                        to_fill += 1
                        pvb = pva + (-1 if p1 == 3 else p1)
                        self.decompression_buffer[to_fill] = pvb
                        if fill_num > 2:
                            to_fill += 1
                            pva = pvb + (-1 if p2 == 3 else p2)
                            self.decompression_buffer[to_fill] = pva
                            if fill_num > 3:
                                to_fill += 1
                                self.decompression_buffer[to_fill] = pva + (-1 if p3 == 3 else p3)
                    self.decompression_prior = self.decompression_buffer[to_fill]
                    dp += fill_num
                    left_to_fill -= fill_num
                cp += 1
            else:
                repeat_count = value - 127
                cp += 1
                repeat_value = self.compression_buffer[self.compression_pos + cp: self.compression_pos + cp+1]
                self.decompression_buffer[self.decompression_pos + dp:
                                          self.decompression_pos + dp + repeat_count] = repeat_value * repeat_count
                dp += repeat_count
                self.decompression_prior = struct.unpack('B', repeat_value)[0]
                cp += 1
        return dp

    def decompress_pbd16(self, compression_len):
        cp = 0
        dp = 0

        def get_pre():
            return struct.unpack(self.endian + 'H',
                                 bytes(self.decompression_buffer[self.decompression_pos + dp - 2:
                                                                 self.decompression_pos + dp]))[0]

        while cp < compression_len:
            code = self.compression_buffer[self.compression_pos + cp]
            if code < 32:
                count = code + 1
                self.decompression_buffer[self.decompression_pos + dp:self.decompression_pos + dp + count * 2] = \
                    self.compression_buffer[self.compression_pos + cp + 1:self.compression_pos + cp + 1 + count * 2]
                cp += count * 2 + 1
                dp += count * 2
                self.decompression_prior = get_pre()
            elif code < 80:
                left_to_fill = code - 31
                while left_to_fill > 0:

                    # 332
                    cp += 1
                    src_char = self.compression_buffer[self.compression_pos + cp]
                    d0 = src_char
                    d0 >>= 5
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', self.decompression_prior + (d0 if d0 < 5 else 4 - d0))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d1 = src_char
                    d1 >>= 2
                    d1 &= 7
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d1 if d1 < 5 else 4 - d1))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d2 = src_char
                    d2 &= 3
                    carry_over = d2

                    # 1331
                    cp += 1
                    src_char = self.compression_buffer[self.compression_pos + cp]
                    d0 = src_char
                    carry_over <<= 1
                    d0 >>= 7
                    d0 |= carry_over
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d0 if d0 < 5 else 4 - d0))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d1 = src_char
                    d1 >>= 4
                    d1 &= 7
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d1 if d1 < 5 else 4 - d1))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d2 = src_char
                    d2 >>= 1
                    d2 &= 7
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d2 if d2 < 5 else 4 - d2))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d3 = src_char
                    d3 &= 1
                    carry_over = d3

                    # 233
                    cp += 1
                    src_char = self.compression_buffer[self.compression_pos + cp]
                    d0 = src_char
                    carry_over <<= 2
                    d0 >>= 6
                    d0 |= carry_over
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d0 if d0 < 5 else 4 - d0))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d1 = src_char
                    d1 >>= 3
                    d1 &= 7
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d1 if d1 < 5 else 4 - d1))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    d2 = src_char
                    d2 &= 7
                    self.decompression_buffer[self.decompression_pos + dp: self.decompression_pos + dp + 2] = \
                        struct.pack(self.endian + 'H', get_pre() + (d2 if d2 < 5 else 4 - d2))
                    dp += 2
                    left_to_fill -= 1
                    if left_to_fill == 0:
                        break
                    self.decompression_prior = get_pre()
                self.decompression_prior = get_pre()
                cp += 1
            elif code < 223:
                raise "unimplemented in vaa3d"
            else:
                repeat_count = code - 222
                cp += 1
                repeat_value = self.compression_buffer[self.compression_pos + cp: self.compression_pos + cp + 2]
                self.decompression_buffer[self.decompression_pos + dp:
                                          self.decompression_pos + dp + repeat_count * 2] = repeat_value * repeat_count
                dp += repeat_count * 2
                cp += 2
                self.decompression_prior = struct.unpack(self.endian + 'H', repeat_value)[0]
        return dp

    def update_compression_buffer8(self):
        look_ahead = self.compression_pos
        while look_ahead < self.total_read_bytes:
            lav = self.compression_buffer[look_ahead]
            if lav < 33:
                if look_ahead + lav + 1 < self.total_read_bytes:
                    look_ahead += lav + 2
                else:
                    break
            elif lav < 128:
                compressed_diff_entries = (lav - 33) // 4 + 1
                if look_ahead + compressed_diff_entries < self.total_read_bytes:
                    look_ahead += compressed_diff_entries + 1
                else:
                    break
            else:
                if look_ahead + 1 < self.total_read_bytes:
                    look_ahead += 2
                else:
                    break
        compression_len = look_ahead - self.compression_pos
        d_length = self.decompress_pbd8(compression_len)
        self.compression_pos = look_ahead
        self.decompression_pos += d_length

    def update_compression_buffer16(self):
        look_ahead = self.compression_pos
        while look_ahead < self.total_read_bytes:
            lav = self.compression_buffer[look_ahead]
            if lav < 32:
                if look_ahead + (lav + 1) * 2 < self.total_read_bytes:
                    look_ahead += (lav + 1) * 2 + 1
                else:
                    break
            elif lav < 80:
                compressed_diff_bytes = int(((lav - 31) * 3 / 8) - 0.0001) + 1
                if look_ahead + compressed_diff_bytes < self.total_read_bytes:
                    look_ahead += compressed_diff_bytes + 1
                else:
                    break
            elif lav < 183:
                compressed_diff_bytes = int(((lav - 79) * 4 / 8) - 0.0001) + 1
                if look_ahead + compressed_diff_bytes < self.total_read_bytes:
                    look_ahead += compressed_diff_bytes + 1
                else:
                    break
            elif lav < 223:
                compressed_diff_bytes = int(((lav - 182) * 5 / 8) - 0.0001) + 1
                if look_ahead + compressed_diff_bytes < self.total_read_bytes:
                    look_ahead += compressed_diff_bytes + 1
                else:
                    break
            else:
                if look_ahead + 2 < self.total_read_bytes:
                    look_ahead += 3
                else:
                    break
        compression_len = look_ahead - self.compression_pos
        d_length = self.decompress_pbd16(compression_len)
        self.compression_pos = look_ahead
        self.decompression_pos += d_length

    def update_compression_buffer3(self):
        pass

    def load_image(self, path):
        assert os.path.exists(path)
        self.decompression_prior = 0
        format_key = "v3d_volume_pkbitdf_encod"
        with open(path, "rb") as f:
            file_size = os.path.getsize(path)
            header_sz = 4 * 4 + 2 + 1 + len(format_key)
            assert (file_size >= header_sz)
            format = f.read(len(format_key)).decode('utf-8')
            assert (format == format_key)
            endian_code_data = f.read(1).decode('utf-8')
            if endian_code_data == 'B':
                self.endian = '>'
            elif endian_code_data == 'L':
                self.endian = '<'
            else:
                raise Exception('endian be big/little')
            datatype = struct.unpack(self.endian + 'h', f.read(2))[0]
            if datatype == 1 or datatype == 33:
                dt = 'u1'
            elif datatype == 2:
                dt = 'u2'
            else:
                raise Exception('datatype be 1/2/4')
            self.load_datatype = datatype
            unit_size = datatype
            sz = struct.unpack(self.endian + 'iiii', f.read(4 * 4))
            self.pbd_sz = sz
            total_unit = sz[0] * sz[1] * sz[2] * sz[3]
            remaining_bytes = compressed_bytes = file_size - header_sz
            self.max_decompression_size = total_unit * datatype
            self.channel_len = sz[0] * sz[1] * sz[2]
            self.pbd3_cur_chan = 0
            read_step_size_bytes = 1024 * 20000
            self.total_read_bytes = 0
            self.compression_buffer = bytearray(compressed_bytes)
            self.decompression_buffer = bytearray(self.max_decompression_size)

            while remaining_bytes > 0:
                current_read_bytes = min(remaining_bytes, read_step_size_bytes)
                self.pbd3_cur_chan = self.total_read_bytes // self.channel_len
                bytes2channel_bound = (self.pbd3_cur_chan + 1) * self.channel_len - self.total_read_bytes
                current_read_bytes = min(current_read_bytes, bytes2channel_bound)
                new_bytes = f.read(current_read_bytes)
                self.compression_buffer[self.total_read_bytes: self.total_read_bytes + current_read_bytes] = new_bytes
                self.total_read_bytes += current_read_bytes
                remaining_bytes -= current_read_bytes
                if datatype == 1:
                    self.update_compression_buffer8()
                elif datatype == 33:
                    self.update_compression_buffer3()
                elif datatype == 2:
                    self.update_compression_buffer16()
            img = np.frombuffer(self.decompression_buffer, self.endian + dt)
            return img.reshape(sz[-1:-5:-1])
