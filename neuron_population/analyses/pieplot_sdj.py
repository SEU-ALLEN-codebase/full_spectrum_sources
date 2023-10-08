ana_types=[th_stypes,cnu_stypes,ctx_stypes]
radar_fea_array=np.zeros((len(ana_types),len(radar_feas)),)
for ity,ty in enumerate(ana_types):
    if ty.count('SSp'):
        ty_df=all_feas[(all_feas.soma_region.isin(ssp_regions))|
                       (all_feas.soma_region.isin(ty))].copy()
    else:
        ty_df=all_feas[all_feas.soma_region.isin(ty)].copy()
    radar_fea_array[ity,0]=ty_df['bnum'].mean(0)
    radar_fea_array[ity,1]=ty_df['tebratio'].mean(0)
    radar_fea_array[ity,2]=ty_df['pd2s_mean'].mean(0)
    radar_fea_array[ity,3]=ty_df['ed2s_mean'].mean(0)
    radar_fea_array[ity,4]=ty_df['binterval_mean'].mean(0)
    radar_fea_array[ity,5]=ty_df['neighborb_mean'].mean(0)
#normalized to 0-100
for i in range(len(radar_feas)):
    radar_fea_array[:,i]/=np.max(radar_fea_array[:,i])
#

radar_fea_df=pd.DataFrame(100*radar_fea_array,columns=radar_feas)
radar_fea_df['B-type']=btype
my_dpi=150
rader_fig=plt.figure(figsize=(4*len(ana_types),4), dpi=my_dpi)
 
# Create a color palette:
my_palette = plt.cm.get_cmap("Set2", len(radar_fea_df.index))
def make_spider(df,row, title,color,type_num=6):
    # number of variable
    categories=list(df)[:-1]
    # print(categories)
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(1,type_num+1,row+1, polar=True, )

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, color='black', size=10)
    # ax.set_xticklabels(rotation=45)

    # Draw ylabels
    ax.set_rlabel_position(0)
    plt.yticks([25,50,75], ["25","50","75"], color="grey", size=9)
    # plt.ylim(0,100)

    # Ind1
    values=df.loc[row].drop('B-type').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, color=color, linewidth=2, linestyle='solid')
    ax.fill(angles, values, color=color, alpha=0.4)
    # Add a title
    plt.title(title, size=16, color='black', y=1.1)
# Loop to plot
for row in range(0, len(radar_fea_df.index)):
    make_spider(radar_fea_df,row=row, title=radar_fea_df['B-type'][row],color=bcolors[row],type_num=len(ana_types))