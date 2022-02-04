# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 21:11:08 2022

@author: njulu
"""

from substrates import SubstrateTab

sub=SubstrateTab()
sub.output_dir='../output'
#sub.plot_svg(1)
sub.plot_substrate(596)





from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio






outputfolder=r"../../PhysiCell-ECM-master/PhysiCell-ECM-master/output"
#outputfolder=r"../../PhysiCell/output"
#outputfolder=r"../../PhysiCell/model1/output"
#outputfolder=r"../../PhysiCell/model2_motility/output"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\test_PhysiCell\PhysiCell_CL\output"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\test_PhysiCell\PhysiCell_drop\output"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\test_PhysiCell\physicellecm-r34\src\output"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\examples_PhysiCell\microenvnmtr-r15\src\output"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\examples_PhysiCell\trmechanics-r7\src\output"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\examples_PhysiCell\trmotility-r17\src"
outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\test_PhysiCell\leader_follower_model\output"
#outputfolder=r"D:\ChengjieLuo\Research_TUE\Friedl Peter and Felix\test_PhysiCell\physicellecm-r34\src\output_1_1"

flag_series=0
tmax=25
dt=1

if flag_series==1:
    first_index=1
    last_index=int(tmax/dt)
    nindex=last_index-first_index+1
    
    
    xs=np.zeros(nindex)
    ys=np.zeros(nindex)
    zs=np.zeros(nindex)
    Nc=np.zeros(nindex)
    
    
    
    
    times=np.zeros(nindex)
    
    f = open(outputfolder+"/cells_series.xyz", "w")
    fd=open(outputfolder+"/cells_series.dump","w")
    for i,n in enumerate(range( first_index,last_index+1,1)):
        try:
            filename='output'+"%08i"%n+'.xml'
            mcds=pyMCDS(filename,outputfolder)
            times[i]= mcds.get_time()
            Nc[i]=mcds.data['discrete_cells']['position_x'].shape[0]
            
            ID=mcds.data['discrete_cells']['ID']
            x=mcds.data['discrete_cells']['position_x']
            y=mcds.data['discrete_cells']['position_y']
            z=mcds.data['discrete_cells']['position_z']
            v=mcds.data['discrete_cells']['total_volume']
            
            boundary=mcds.data['mesh']['boundary']
            xmin=boundary[0]
            xmax=boundary[3]
            ymin=boundary[1]
            ymax=boundary[4]
            zmin=boundary[2]
            zmax=boundary[5]
            
            # fig = plt.figure()
            # ax = fig.add_subplot(projection='3d')
            # ax.scatter(x,y,z)
            
            # plt.xlim([-100,100])
            # plt.ylim([-100,100])
            # plt.xlim([min(x)-1,max(x)+1])
            # plt.ylim([min(y)-1,max(y)+1])
            
            # plt.figure()
            # plt.plot(x,y,'o')
          
            # # cycle=mcds.data['discrete_cells']['cycle_model']
            # # p = mcds.data['discrete_cells']['oncoprotein']
            # x=mcds.data['discrete_cells']['position_x']
            # y=mcds.data['discrete_cells']['position_y']
            # z=mcds.data['discrete_cells']['position_z']
        
            # xs[i]=x
            # ys[i]=y
            # zs[i]=z
            
            ## write to xyz format
            f.write(str(int(Nc[i]))+'\n')
            f.write('Atoms. Timestep:'+str(int(times[i]))+'\n')
            for ic in np.arange(int(Nc[i])):
                ###### write to xyz
                tmpx=str(x[ic])
                tmpy=str(y[ic])
                tmpz=str(z[ic])
                tmpid=str(ID[ic])
                tmpv=v[ic]
                tmpr=str((tmpv/(4/3.*np.pi))**(1/3))
                
                f.write(tmpid+' '+tmpx+' '+tmpy+' '+tmpz+' '+tmpr+'\n')
                
                
            ###### write to dump format
            fd.write("ITEM: TIMESTEP \n")
            fd.write(str(int(times[i]))+'\n')
            fd.write("ITEM: NUMBER OF ATOMS \n")
            fd.write(str(int(Nc[i]))+'\n')
            fd.write("ITEM: BOX BOUNDS pp pp pp\n")
            fd.write(str(xmin)+' '+str(xmax)+ "\n")
            fd.write(str(ymin)+' '+str(ymax)+ "\n")
            fd.write(str(zmin)+' '+str(zmax)+ "\n")
            fd.write("ITEM: ATOMS type id x y z radius \n")
            for ic in np.arange(int(Nc[i])):
                ###### write to xyz
                tmpx=str(x[ic])
                tmpy=str(y[ic])
                tmpz=str(z[ic])
                tmpid=str(int(ID[ic]))
                tmpv=v[ic]
                tmpr=str((tmpv/(4/3.*np.pi))**(1/3))
                
                fd.write(tmpid+' '+'1'+' '+tmpx+' '+tmpy+' '+tmpz+' '+tmpr+'\n')
        except:
            print("not all included")     
            break
    
    
    
    
    f.close()
    fd.close()
    
    print(Nc)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x,y,z)
    plt.show()


n=200
filename='output'+"%08i"%n+'.xml'
#filename='initial.xml'
mcds=pyMCDS(filename,outputfolder)




ID=mcds.data['discrete_cells']['ID']
x=mcds.data['discrete_cells']['position_x']
y=mcds.data['discrete_cells']['position_y']
z=mcds.data['discrete_cells']['position_z']
v=mcds.data['discrete_cells']['total_volume']
cycle_model=mcds.data['discrete_cells']['cycle_model']
cell_type=mcds.data['discrete_cells']['cell_type']
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(x,y,z)

fig = plt.figure()
plt.plot(x,y,'o')


leaderx=x[cell_type==1]
leadery=y[cell_type==1]
followerx=x[cell_type==2]
followery=y[cell_type==2]
fig = plt.figure()
plt.plot(leaderx,leadery,'bo',markeredgecolor='k')
plt.plot(followerx,followery,'yo',markeredgecolor='k')


mcds.get_substrate_names();
#o2 = mcds.get_concentrations( 'substrate' );
o2 = mcds.get_concentrations( 'oxygen' );
#o2 = mcds.get_concentrations( 'ECM anisotropy' );

X,Y = mcds.get_2D_mesh();
plt.figure()
plt.contourf(X,Y,o2[:,:,0],20);
#plt.plot(x,y,'r.')
#plt.contourf(X,Y,o2[:,:,10]);

plt.colorbar()
plt.axis('image')




################ plot ECM ##################

filename=outputfolder+'/output'+"%08i"%n+'_ECM.mat'
ecm=sio.loadmat(filename)['ECM_Data']
centers=ecm[:2,:]
ecmanisotropy=ecm[3,:].reshape(75,75)
fiber_alignment=ecm[5:7,:]



plt.figure()
plt.contourf(X,Y,ecmanisotropy,20,cmap='YlOrRd');
#plt.plot(x,y,'r.')
#plt.contourf(X,Y,o2[:,:,10]);

plt.colorbar()
plt.axis('image')


#plt.figure()
#plt.plot(centers[0,:],centers[1,:],'.')
xml_fname = "output%08d.xml" % n
snapshot = xml_fname[:-4]
mcds.load_ecm(snapshot + '_ECM.mat', outputfolder)

cell_df = mcds.get_cell_df()
xx, yy = mcds.get_2D_mesh()
micro = mcds.get_concentrations('ECM anisotropy', 0.0)

# find levels for microenvironment
# plane_oxy = mcds.get_concentrations('oxygen', 0.0)
# num_levels = 25
# #levels = np.linspace(plane_oxy.min()+1e-14, plane_oxy.max(), num_levels)
# levels = np.linspace(1e-14, 38, num_levels)

# arrow lengths depend on anisotropy
micro_scaled = micro
#print_stats(micro_scaled)
#mean = np.mean(micro_scaled.flatten())
V_max = 4
#K_M = mean
K_M = 0.4
def curve(x):
    #return (V_max * x) / (K_M + x)
    return 0.5 if x > 0.5 else x

for i in range(len(micro)):
    for j in range(len(micro[i])):
        #micro_scaled[i][j] = 10 *  math.log10(micro[i][j] + 1) / math.log10(2)
        micro_scaled[i][j] = curve(micro[i][j])
        
micro_scaled = micro
#print_stats(micro_scaled)

dy = mcds.data['ecm']['y_vec'][:, :, 0] * micro_scaled
dx = mcds.data['ecm']['x_vec'][:, :, 0] * micro_scaled
#print(dx.shape)
#print('dmag (min, max)', (np.sqrt(dx**2 + dy**2).min(), np.sqrt(dx**2 + dy**2).max()))

# normalize lengths -- this needs some fiddling
#dx = dx / dx.std()
#dy = dy / dy.std()

# if we want the arrows the same length instead
dx_unscaled = mcds.data['ecm']['x_vec'][:, :, 0]
dy_unscaled = mcds.data['ecm']['y_vec'][:, :, 0]

# mask out zero vectors
mask = np.logical_or(dx > 1e-4, dy > 1e-4)

# add quiver layer with scaled arrows ###
plt.quiver(xx[mask], yy[mask], dx[mask], dy[mask], pivot='middle', angles='xy', units='width', headwidth=0, width=.0015)
plt.xlim([-400,400])
plt.ylim([-400,400])
#plt.quiver(xx[mask], yy[mask], dx_unscaled[mask], dy_unscaled[mask], pivot='middle', angles='xy', units='width', headwidth=0, width=.0015)







# x10=-500
# y10=500
# x1=x[(x<0)&(y>0)]
# y1=y[(x<0)&(y>0)]
# x10=x1.mean()
# y10=y1.mean()
# r1=np.sqrt((x1-x10)**2+(y1-y10)**2)
# plt.scatter(x1,y1)

# x20=-500
# y20=-500
# x2=x[(x<0)&(y<0)]
# y2=y[(x<0)&(y<0)]
# x20=x2.mean()
# y20=y2.mean()
# r2=np.sqrt((x2-x20)**2+(y2-y20)**2)
# plt.scatter(x2,y2)

# x30=500
# y30=500
# x3=x[(x>0)&(y>0)]
# y3=y[(x>0)&(y>0)]
# x30=x3.mean()
# y30=y3.mean()
# r3=np.sqrt((x3-x30)**2+(y3-y30)**2)
# plt.scatter(x3,y3)


# x4=x[(x>0)&(y<0)]
# y4=y[(x>0)&(y<0)]
# x40=x4.mean()
# y40=y4.mean()
# r4=np.sqrt((x4-x40)**2+(y4-y40)**2)
# plt.scatter(x4,y4)




# plt.figure()
# for r in [r1,r2]:
#     a1,b1=np.histogram(r,bins=30,range=[0,500],density=True)
#     plt.plot((b1[:-1]+b1[1:])/2/r.max(),a1,'o-',label='slow')
# for r in [r3,r4]:
#     a1,b1=np.histogram(r,bins=30,range=[0,500],density=True)
#     plt.plot((b1[:-1]+b1[1:])/2/r.max(),a1,'s--',label='fast')
    
# plt.xlabel('normed r')
# plt.ylabel('PDF')
# plt.legend()


# plt.figure()
# for r in [r1,r2]:
#     a1,b1=np.histogram(r,bins=30)
#     plt.plot((b1[:-1]+b1[1:])/2,a1,'o-',label='slow')
# for r in [r3,r4]:
#     a1,b1=np.histogram(r,bins=30)
#     plt.plot((b1[:-1]+b1[1:])/2,a1,'s--',label='fast')
    
# plt.xlabel('r')
# plt.ylabel('N(r)')
# plt.legend()
