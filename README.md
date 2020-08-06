# FLUID SIMULAITION ON PROCESSING + AUDIO VISUALIZER
Inspired by Daniel Shiffman's work about Fluid Simulation, coding challenge #132. Thanks a lot for your amazing content, sincerly.

# INTRO
I've always been amazed by fluids' behaviour, and it's been a while since I thought to try to write some code about it.
My goal is to learn shading languages and implement the algorithm on GPU in order to unleash the full potential of this kind of simulation, but it's not my cup of tea at the moment so CPU is the way to go. I even bought GPU Gems (https://www.nvidia.com/en-us/drivers/gpu-gems-home/) in order to understand something more and have a good guide to follow, I just need to put it in works really.

# DEPENDECIES
In order to run the simulation you must install Processing (last version is 3.5.4, currently). Then, inside Processing, you have to install two libraries through Sketch->Import Library->Add Library:
- ControlP5
- Sound

# HOW TO USE
If you don't want to use the audio visualization, just remove "audio" from methods in draw() (look in fluid_simulation.pde); otherwise, you must create a folder in the main directory called "data" and put in this folder the track you want to visualize, in .wav format. Then, you can change the string "bias.wav" in setup() (fluid_simulation.pde) with the exact name of the track you want to visualize.

During simulation, if you're using it in "audio" mode, you can change density and velocity multipliers value by tweaking "dMul" and "vMul".
