#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : feature_processing.py
#   Author       : Yufeng Liu
#   Date         : 2022-09-01
#   Description  : 
#
#================================================================

def whitening(features, epsilon=1e-9):
    # do whitening
    features = (features - features.mean(axis=0)) / (features.std(axis=0) + epsilon)
    return features
