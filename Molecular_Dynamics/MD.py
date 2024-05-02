import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#-----------------------------Question 2 & 3----------------------------------------------
#question2: (time, KE, PE, Tot, 20K iter, data 10 itr, no thermo, no eqblm, dt=0.005)
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q2_energies_dt(0.005).dat') 

#question3: (time, KE, PE, Tot, 20K iter, data 10 itr, thermo at 1000 itr, no eqblm, dt=0.005)
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q3_energies_TS(1000).dat') 

# N=2197 #no. of particles
# x = data[:,0]*0.005 #time
# PE = data[:,1]
# KE = data[:,2]
# TOT = data[:,3]
# print('For dt = 0.005')
# print('PE mean and std dev', round(np.mean(PE[800:]),4), round(np.std(PE[800:]),4))
# print('KE mean and std dev', round(np.mean(KE[800:]),4), round(np.std(KE[800:]),4))
# print('E mean and std dev', round(np.mean(TOT[800:]),4), round(np.std(TOT[800:]),6))

# plt.plot(x,PE,'-',markersize=0.5)
# plt.plot(x,PE,'.',markersize=4,label='PE')

# plt.plot(x,KE,'-',markersize=0.5)
# plt.plot(x,KE,'.',markersize=4,label='KE')

# plt.plot(x,TOT,'-',markersize=0.5)
# plt.plot(x,TOT,'.',color='k',markersize=4,label='Tot')

# #change as per question:
# plt.suptitle('Energy variation with time')
# plt.title('2197 particles, dt=0.005')
# plt.ylabel('Energy')
# plt.xlabel('Time')
# plt.legend()
# # plt.ylim(0.30755,0.30762)
# # plt.savefig('q3_energy_dt(0.005).png')
# plt.show() 


# #momentum conservation
# data1 = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q3_velocity_dt(0.0025).dat') 
# x = data[:,0]*0.0025 #time
# vx = data1[:,0]
# vy = data1[:,1]
# vz = data1[:,2]

# plt.plot(x,vx,'.',markersize=4,label='Px')
# plt.plot(x,vy,'.',markersize=4,label='Py')
# plt.plot(x,vz,'.',color='k',markersize=4,label='Pz')

# plt.suptitle('Momentum variation with time')
# plt.title('2197 particles, dt=0.0025')
# plt.ylabel('Momentum')
# plt.xlabel('Time')
# plt.legend()
# # plt.ylim(0.30755,0.30762)
# # plt.savefig('q3_momentum_dt(0.0025).png')
# plt.show() 


#-----------------------------Question 4 & 5----------------------------------------------
#question 4: neigh_list:
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q4_neighbour_list.dat')

#         # histogram
# plt.hist(data, bins=20, edgecolor='black', color='skyblue')
# plt.xlabel('Number of neighbours')
# plt.ylabel('Frequency')
# plt.suptitle('Most number of neighbours in every iteration',fontsize=14)
# plt.title('Between 20K-30K iterations, 1200 particles',fontsize=9)
# # plt.savefig('q4_hist_neighlist.png')
# plt.show()

# #question 5: g(r)
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q5_gr.dat') 
# r = data[:,0] #Interparticle distance
# gr = data[:,1]
# print(max(gr))
# plt.plot(r,gr,'.-')
# plt.axvline(x=(max(r)/2),linestyle='--', color='red', label='$r_{max}$')
# plt.xlabel('Interparticle distance ($r$)')
# plt.ylabel('$g(r)$')
# plt.suptitle('Pair Correlation function',fontsize=14)
# plt.title('1200 particles, dr=0.1$\sigma$',fontsize=9)
# plt.text(2, 2, f'Peak height:{round(max(gr),3)}', fontsize=12, color='k')
# plt.xlim(0,11)
# plt.legend()
# plt.grid(True)
# # plt.savefig('q5_gr.png')
# print(f'Peak height:{round(max(gr),3)}')
# plt.show()

# #-----------------------------Question 6----------------------------------------------
# # n_part in Q6=2400
# #neigh_list:
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q6_neighbour_list.dat')

#         # histogram
# plt.hist(data, bins=20, edgecolor='black', color='skyblue')
# plt.xlabel('Number of neighbours')
# plt.ylabel('Frequency')
# plt.suptitle('Most number of neighbours in every iteration',fontsize=14)
# plt.title('Between 20K-30K iterations, 2400 particles',fontsize=9)
# # plt.grid(True)
# plt.show()

# #g(r)
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q6_gr.dat') 
# r = data[:,0] #Interparticle distance
# gr = data[:,1]
# print(max(gr))

# new_r = []
# new_gr =[]
# for i in range(len(r)):
#     if data[i,0] > 2 and data[i,0] <3:
#         new_r.append(data[i,0])
#         new_gr.append(data[i,1])

# print([new_r])
# print([new_gr])

# plt.plot(r,gr,'.-')
# plt.axvline(x=(max(r)/2),linestyle='--', color='red', label='$r_{max}$',ymax=0.5)
# # plt.axvline(x=(new_r[2]),linestyle='--', color='k', label='$r_{max}$',ymax=0.6,ymin=0)
# plt.xlabel('Interparticle distance ($r$)')
# plt.ylabel('$g(r)$')
# plt.suptitle('Pair Correlation function',fontsize=14)
# plt.title('2400 particles, dr=0.1$\sigma$',fontsize=9)
# # plt.arrow(2,2, 1, round(max(gr),3), head_width=0.1, head_length=0.2, fc='red', ec='black')
# plt.text(2,2, f'First Peak height:{round(max(gr),3)}', fontsize=12, color='k')
# plt.text(3, 1.8, f'Second Peak height:{round(max(new_gr),3)}', fontsize=12, color='k')
# # plt.text(4, 1.6, f'Third Peak height:{round(sorted_y_values[2],3)}', fontsize=12, color='k')
# plt.xlim(0,11)
# plt.grid(True)
# plt.legend()
# print(f'Second Peak height:{round(max(new_gr),4)}')
# plt.savefig('q6_gr.png')
# plt.show()

# #-----------------------------Question 7----------------------------------------------
# # n_part in Q7=3600
# #neigh_list:
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q7_neighbour_list.dat')

#         # histogram
# plt.hist(data, bins=20, edgecolor='black', color='skyblue')
# plt.xlabel('Number of neighbours')
# plt.ylabel('Frequency')
# plt.suptitle('Most number of neighbours in every iteration',fontsize=14)
# plt.title('Between 20K-30K iterations, 3600 particles',fontsize=9)
# # plt.grid(True)
# plt.show()

# #g(r)
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q7_gr.dat') 
# r = data[:,0] #Interparticle distance
# gr = data[:,1]
# print(max(gr))

# new_r2 = []
# new_gr2 =[]
# new_gr3 =[]
# new_r3= []
# for i in range(len(r)):
#     if data[i,0] > 2 and data[i,0] <3:
#         new_r2.append(data[i,0])
#         new_gr2.append(data[i,1])
#     if data[i,0] > 3 and data[i,0] <4:
#         new_r3.append(data[i,0])
#         new_gr3.append(data[i,1])

# print([r[:15]])
# print([gr[:30]])
# print([new_r2])
# print([new_gr2])
# print([new_r3])
# print([new_gr3])

# plt.plot(r,gr,'.-')
# plt.axvline(x=(max(r)/2),linestyle='--', color='red', label='$r_{max}$',ymax=0.5)
# # plt.axvline(x=(new_r[2]),linestyle='--', color='k', label='$r_{max}$',ymax=0.6,ymin=0)
# plt.xlabel('Interparticle distance')
# plt.ylabel('$g(r)$')
# plt.suptitle('Pair Correlation function',fontsize=14)
# plt.title('3600 particles, dr=0.1$\sigma$',fontsize=9)
# plt.text(2, 1.8, f'First Peak height:{round(max(gr),3)}', fontsize=12, color='k')
# plt.text(3, 1.6, f'Second Peak height:{round(max(new_gr2),3)}', fontsize=12, color='k')
# plt.text(4, 1.4, f'Third Peak height:{round(max(new_gr3),3)}', fontsize=12, color='k')
# plt.xlim(0,11)
# plt.grid(True)
# plt.legend()
# plt.savefig('q7_gr.png')
# print(f'First Peak height:{round(max(gr),4)}')
# print(f'Second Peak height:{round(max(new_gr2),4)}')
# print(f'Third Peak height:{round(max(new_gr3),4)}')
# plt.show()

# #-----------------------------Question 8----------------------------------------------
# import math
# data1 = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\q8_vel_dist_3600p.dat') 
# # data2 = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\Dat files\q8_vel_dist_2400p.dat') 
# # data3 = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\Dat files\q8_vel_dist_3600p.dat') 



# plt.plot(data1,'.-',markersize=7,linewidth=1)#, label='1200')
# # plt.plot(data2,'.-',markersize=7,linewidth=1,label='2400')
# # plt.plot(data3,'.-',markersize=7,linewidth=1,label='3600')
# plt.xlabel('Velocity ($|v_x|$)')
# plt.ylabel('$P(v)$')
# plt.suptitle('Velocity distribution',fontsize=14)
# # plt.title('For three set of particles',fontsize=9)
# plt.xlim(0,100)
# plt.legend()
# plt.grid(True)
# # plt.savefig('q8_MBD.png')
# plt.show()


#----------------------------- Particle in 3D box
# num_particles = 600
# L = 20 # Replace L with box size
# data = np.loadtxt('G:\Fortran_workspace_CompPhy\Assignment7_MD\initial_position.dat',skiprows=1)
# x_positions = data[:,0]  
# y_positions = data[:,1]
# z_positions = data[:,2]

# # Create a 3D scatter plot
# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x_positions, y_positions, z_positions, s=10, color='b', marker='o')

# ax.set_xlim(0, L)
# ax.set_ylim(0, L)
# ax.set_zlim(0, L)
# ax.set_xlabel('X position')
# ax.set_ylabel('Y position')
# ax.set_zlabel('Z position')
# ax.set_title('Particle Positions in a 3D Box')
# plt.show()
