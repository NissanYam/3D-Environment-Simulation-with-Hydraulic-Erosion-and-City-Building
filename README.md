# 3D Environment Simulation with Hydraulic Erosion and City Building

This project is a simulation of a 3D environment where the user can navigate through a terrain, simulate hydraulic erosion, and build a city. The simulation is developed using C++ and OpenGL, utilizing the GLUT library for window management and input handling.

## Prerequisites

- C++ compiler that supports C++11 (e.g., g++)
- OpenGL
- GLUT library

## How to Build

1. Compile the source files using a C++ compiler:

    ```bash
    g++ -o main src/main.cpp -lglut -lGLU -lGL -lm
    ```

2. Run the compiled executable:

    ```bash
    ./main
    ```

## Controls

- **Arrow Keys**: Change the direction of movement.
- **Page Up/Page Down**: Move the camera up/down.
- **Left/Right Arrow Keys**: Rotate the camera left/right.
- **Up/Down Arrow Keys**: Increase/decrease movement speed.
- **Mouse Left Click**: Start/stop rain simulation.

## Menu Options

1. **Regular view**: Switches to the regular view mode.
2. **Top view**: Switches to the top-down view mode.
3. **Build city**: Initiates the city building process.
4. **Reset world**: Resets the world to its initial state.
5. **Reset to starting world**: Resets the world to its starting state.

## Features

1. **Navigation**: User can move around the environment and look in different directions.
2. **Rain Simulation**: Rainfall can be simulated, causing hydraulic erosion.
3. **City Building**: A city can be built on the terrain.

## How it Works

1. The program creates a 3D terrain with randomly generated elevations.
2. The user can navigate through the terrain using keyboard controls.
3. Rainfall can be simulated, eroding the terrain and creating valleys.
4. City building involves selecting points on the terrain and using a flood-fill algorithm to construct buildings.
5. Menu options allow for switching between views, starting/stopping rain, and building/resetting the city.

## Credits

This project was created by Nissan Yamin.
