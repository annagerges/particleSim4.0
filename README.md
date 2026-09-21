# Unified Multi-Physics Particle System

A C++ simulation modeling gravitational and spring forces on particles in 2D space, with spatial partitioning for efficient collision detection and RK4 numerical integration for physics accuracy. Python for data analysis.

---

## Overview

This project simulates a system of 100–1000 particles affected by gravity and spring forces, with particle-to-particle and wall collisions. The simulation uses:

- **RK4 integration** for <0.1% energy drift and increase in accuracy compared to Euler integration
- **Spatial partitioning hashmap** to reduce collision detection from O(n²) to O(n)
- **CSV export** for post-simulation analysis and validation
- **Python validation plots** to visualize energy conservation

---

## Physics Model

### Forces and Acceleration

**Gravitational Force (freefall)**:

```
a = -g = -9.8 m/s²
```

**Spring Force (when y ≤ spring height)**:
I wanted to  experiment with repulsive forces to keep the simulation lively.

```
F_spring = k(y - h)
a_total = (k/m)(y - h) - g
acceleration is the force of the spring/mass minus gravity

The spring stiffness `k` is scalable as a cmd line parameter but if no arguments are provided it is 
dynamically calculated to support a particle at equilibrium compression (0.2 m) where m is the mass per particle:
k = ( 9.8 × m / 0.2) × 4

Factor of 4 increases springiness and creates a responsive force.
```

### Numerical Integration

**RK4 (Runge-Kutta 4th Order)** with fixed timestep `dt = 0.1 s`. Computes four slope estimates (k1–k4) and combines them with weights (1:2:2:1) to approximate the solution:
Though RK4 take more computational power that Euler, it is magnitudes more accurate over time.

```cpp
		// Determine if particle is on spring once at start of step
		bool onSpring = (p.getY() <= s.getHeight());

		// RK4 stages
		float k1y = p.getVy();

		//if it's on the spring change it's accelaration to k*dy/m-g and if not than keep accelaration to -9.8
		float k1v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - p.getY()) - 9.8f : -9.8f;

		float k2y = p.getVy() + 0.5f * k1v * dt;
		float k2v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + 0.5f * k1y * dt)) - 9.8f : -9.8f;

		float k3y = p.getVy() + 0.5f * k2v * dt;
		float k3v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + 0.5f * k2y * dt)) - 9.8f : -9.8f;

		float k4y = p.getVy() + k3v * dt;
		float k4v = onSpring ? (s.getK() / p.getMass()) * (s.getHeight() - (p.getY() + k3y * dt)) - 9.8f : -9.8f;

		// Update to a weighted avg of every slope
		p.setY(p.getY() + (dt / 6.0f) * (k1y + 2 * k2y + 2 * k3y + k4y));
		p.setVy(p.getVy() + (dt / 6.0f) * (k1v + 2 * k2v + 2 * k3v + k4v));

```
Updates apply at (`0.01 s` intervals) within each frame for stability. Rk4 achieves <0.1% energy drift over long simulations compared to the approximately 2% drift with Euler integration.

## Program Architecture

### Spatial Partitioning

The simulation uses a dynamically sized to reduce collision detection from O(n²) to O(n) for typical particle distributions:
-**Grid Sizing**: Width and height are the maximum between 3 and  grid(sqrt(num of particles/10)+2)
- **Grid cell width**: 800 / cells per row
- **Collision checks**: Only particles in the same grid cell are tested
- **Grid update**: `fix()` partial rebuild if any particles move cells to update the simulation while managing resources
- **Multiplier**: Creates a unique key using (row*cellsPerGrid)+cell

This optimization scales efficiently to 1000 particles without performance degradation.

### Collision Detection

#### Particle–Particle Collisions
- **Trigger**: Euclidean distance `absDx * absDx + absDy * absDy < 100` units (r^2) (are the particles overlapping?)
- **Axis determination**: Compare `|dx|` vs. `|dy|` to resolve collision along correct axis
	- **Response**:
	  - **Y-axis collision**: Swap vy for energy conservation
	  - **X-axis collision**: Swap vx for energy conservation
	  - **If X and Y axis are equal**: apply Y-axis collision conditions.

#### Wall Collisions
- **Bounds**: `y values are from 6-max height and x values are from 0-max height (non-inclusive)`
- **Response**: Reverse normal velocity component, clamp position to boundary

## CSV Logging
Uses fstream C++ library to create a csv file and log the qualities of every particle (besides mechanical energy and mass (automatically 0.5kg to prevent division by 0 errors) for further analysis.

#### Python CSV Analysis
- **Method**: Reads CSV, computes kinetic, gravitational potential, and potential spring energy. Plots energy drift % over time. Also let's user choose how many particle trajectories they want to view.
- **Result**: RK4 method has <0.1% energy drift and demonstrates numerical stability
- **Energy Conservation**: ![Energy Conservation Demo](EnergyConservation.png)
- **Particle Trajectories**: ![Particle Trajectory Demo](Trajectory.png)
  
  

## CSV Error Handling
Safely processes invalid data during analysis

```
##go through every row
for _, row in df.iterrows():
    #extract particle number
    raw_pid=row['Particle num']

    try:
        ##check if particle num is valid and convert it to an integer
        p_id=int(float(raw_pid))
    except (ValueError, TypeError):
        p_id=None

    ##if the particle number is a number and is one of the target particles convert the number to an integer, the coordinates to number (with NaN error handling), and append it to the dictionary if the coordinates are numbers
    if p_id is not None and p_id in target_particles:
        px=pd.to_numeric(row['x'], errors='coerce')
        py=pd.to_numeric(row['y'], errors='coerce')

        if not (pd.isna(px) or pd.isna(py)):
            particle_paths[p_id]['x'].append(px)
            particle_paths[p_id]['y'].append(py)
```
Prevents crashes and excludes invalid data.
	

### Frame Timing

Uses `std::chrono::high_resolution_clock` to decouple frame rate from simulation timestep:
```
accumulator += frameDuration.count()
while (accumulator >= dt) {
    update()
    accumulator -= dt
}
```

Ensures consistent physics independent of frame rate or system load.

## Code Structure

| File | Purpose |
|------|---------|
| `particleSim2.0.cpp` | Initialization, user input validation, grid setup, main loop |
| `Particles.cpp` | Physics update, collision detection, grid management |
| `Particles.h` | Class definitions (`Particles`, `Spring`), function declarations |

### Key Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `updatePos()` | `void(vector<Particles>&, Spring&)` | Apply forces, update velocity and position |
| `wallCollis()` | `void(vector<Particles>&)` | Handle boundary collisions |
| `particleCollis()` | `void(std::unordered_map<int, std::vector<Particles*>>&)` | Detect and resolve particle–particle collisions within grid cells with more than 1 particle |
| `fix()` | `void(std::vector<Particles>&, std::unordered_map<int, std::vector<Particles*>>&, int, int)` | Clears and puts any particle that moves cells in accordance to it's current position while maintaining O(1) time|
| `csvDump()` | `void csvDump(std::vector<Particles>&, std::fstream&);` | Logs particle state into csv file for futher analysis|

### Spring Stiffness Tuning
In cmd line the user can set k to a desired number if they would like.

If user enters too many arguments, k will be dynamically calculated instead as to not break the simulation entirely. If the argument is too small or not enough arguments were entered, the program will let the user know and terminate.
```
    if (argc > 2) {
        //set k assuming all of the particles are statically laying on the spring and compressing 0.2m. Mulitply by 4 so that the particles are springy.
        s.setK(((nP * 9.8 * particles[0].getMass()) / 0.2) * 4);
        cout << "Too many arguments entered, so k was dynamically allocated to: " << s.getK();
    }
    else if (argc == 2) {
        double k = stof(argv[1]);
        if (k <=0) {
            cout << "INVALID ARGUMENT. K is too small.";
            return 1;
        }
        else {
            //customize k
            s.setK(k);
        }


        cout << "Argument detected! k dynamically allocated to: " << s.getK() << "\n";
    }
    else {
        // If they clicked the normal VS button (argc == 1) or passed too many arguments (argc > 2)
        cout << "Error: You must provide a command-line argument to run this simulation.\n";
        return 1;
    }
```

#### Requirements
- **C++**: C++11 or later
- **Python**: 3.14+, with `pandas` and `matplotlib`


## Coming Soon!!
- **Rendering**
- **Performance benchmarking**: Analyzing time complexity
- **Damping and Spring Stiffness CMD line Parameters**

---

**Date**: August 2026  
**License**: MIT
