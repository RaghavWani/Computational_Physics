## Molecular Dynamics
The main code is named 'MD.f90', while other .f90 files are subroutine files used in the main code. 
### Procedure
1. Calculate the initial position and velocity.
2. Calculate forces between particles.
3. Update position and velocity using the *Verlet* algorithm.
4. Calculate updated position and velocity.
5. Continue this process over several iterations to get the evolution of the system.

### Subroutines
- **post_init:** It is a subroutine to initialize the position of a specified number of particles. Initialization is done by randomly assigning coordinates to each particle in a 3D box, considering particles are uniformly distributed throughout the box. Here, we ensure that the particles are at least one sigma far (particle diameter) so that no two particles overlap. We store the position data of all particles in a single column such that each set of 3 rows corresponds to individual particles: the first line is x coord, the second line is y coord, and the third line is z coord.

- **vel_init:** We consider $v=A(r - 0.5)$ as the velocity of an individual particle, where $r$ is a random number between 0 and 1. Considering the equipartition of energy per degree of freedom, that is $0.5K_bT/m=0.5 \left\langle{V_x^2}\right\rangle$. After solving the integral, we get $A=\sqrt{12} (K_bT/m)$. substituting this value in v, we get the velocity of each particle along three directions. Next, we do overall velocity correction by subtracting the velocity of the system's center of mass from the velocity of each particle in each direction.

- **calc_force:** We consider Lennard-Jones potential to get force. To make the problem computationally solvable, we modify the potential such that both force and potential go to zero at some $r=r_c$ and beyond. Thus, while calculating the force, we only consider particles within radius $r_c$. Hence, we get modified force and potential as:
$$F'(r) = F(r) - F_c$$
$$V''(r) = V(r) + F_c r - V_c \text{ where, }$$

$$F(r) = 4\epsilon [ 12 \frac{\sigma^6}{r^{13}} - 6 \frac{\sigma^6}{r^{7}}]$$
$$V(r) = 4\epsilon [ (\frac{\sigma}{r})^{12} - (\frac{\sigma}{r})^{6}]$$
$$F_c = F(r) | _{r=r_c} \text{ ; } V_c = (V(r) + F_c r)| _{r=r_c}$$ 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Here, $\sigma$ is particle diameter and $\epsilon$ is potential depth. We have also used the minimum image convention, which corrects the distance between two particles considering periodic boundary conditions. So, only if the distance $x$ between particles is greater than half of the box size ($L/2$), we correct the distance to$$L - x$, or else we consider the distance to be $x$.

- **update_post:** We use the *Verlet* algorithm to update the position of each particle. Formula:
$$x(t+dt) = x(t) + v(t)dt + 0.5 a(t)dt^2$$

- **update_vel:** Similar to position, we update the velocity. But here, we also need force at the current time $t \text{+}dt$ and the previous time $t$. Hence, we calculate the force using new position coordinates of particles and use this force to get velocity at $t \text{+}dt$:
$$v(t+dt) = v(t) + 0.5(\frac{F(t) + F(t+dt)}{2m})$$

- **thermostat:** We observe that over iterations, the system's energy increases and deviates from the classically expected value (i.e., $0.5K_B T$) due to the increasing kinetic energy of particles. So, we rescale kinetic energy to the appropriate value at fixed intervals. Thus, the thermostat takes some energy from the system and lowers its overall energy. This is done until equilibrium is reached. 

The entire process runs over several iterations, and physical quantities like kinetic energy, potential energy, and total energy are calculated regularly. Plotting it against time helps to study the dynamics of the system.

A faster method to do this process is using a neighboring list. While calculating force, we consider the contribution of all the other particles. However, the far enough particles may not significantly affect the force. Hence, we only consider particles that lie within $r=r_s$ around each particle and list them after periodic intervals.

- **neigh_list_finder:** This subroutine is used to get the neighboring list, a set of particles within some fixed radius $r_s$ from each particle. Hence, each particle has a unique neighbor list, which is updated after every desired number of iterations.

- **calc_force_neighb:** It is the same as the *calc_force* subroutine; it calculates contributions only from neighbor list particles.

- **calc_gr:** This subroutine calculates the pair correlation function $g(r)$ between particles. The pair correlation function is defined as the probability that a pair of particles are separated by a distance, say $r$. The plot of $g(r)$ vs. $r$ can help identify the density of particles within the box and distinguish the physical phase of the system (higher-density systems have longer oscillatory behavior)

- **speed_dist:** This subroutine calculates the velocity of particles and gets the data of the speed of particles between $v$ and $v \text{+}dv$. Plotting this data gives us Maxwell-Boltzmann distribution, as expected from ideal particles. One can also plot velocity distribution, which would be a Gaussian.

###### *Summarized from: https://www.youtube.com/playlist?list=PLyqSpQzTE6M8Lg4pPC2KKutByiFCR9kV0*

