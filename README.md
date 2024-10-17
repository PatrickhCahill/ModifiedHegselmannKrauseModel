# ModifiedHegselmannKrauseModel
This repository is designed to reproduce the simulated figures from the *A modified Hegselmann–Krause model for interacting voters and political parties* (2024) by Patrick Honnery Cahill and Georg Gottwald.

## Minimal Reproduction
To reproduce the figures [juliaLANG needs to be installed](https://julialang.org/downloads/) - bash will also be required. Now, open a terminal and execute the following:
* `mkdir data` This will create a data folder where the raw data is stored. Automatically gitignored.
* `julia requirements.jl` This will install the necessary packages in order to execute the code.
* `julia 4to8pipeline.jl` This populate the figures from 4 to 8.
* `julia voter_behaviours_explanation.jl` This will create each subfigure in Figure 2 from the paper.

Now open `figure3.jl` in an editor and edit the foldername variable to correspond to the folder of the "party_base" simulation - which will typically of the form `{todays_data}_run_8` if the above files are executed. With this variable assigned run the following:
* `julia figure3.jl` - This will create figure 3.
* `phasediagrams.bash` will create execute the heatmap related julia files in parallel in order to create figures 10 and 11. **Note requires bash and is computationally expensive.**


Note: Figures are provided with a seed to be reproducible. But they may differ between versions of julia and package versions. The figures in text were produced using `julia version 1.11.0-beta1`.
## sim/Dorder_run
This folder with corresponding data is the output of `phasediagrams.bash`. This corresponds to figures 10 and 11 in text.
## Code explainer
For every simulation we require the instantiation of set of parameters called `the_params` which is an instance of the `Params` struct. Here, we define the interaction forces, interaction radii, total time, time step size, boundary_conditions etc.

We then handle the creation of the initial conditions using the `generate_initial_conditions(the_params, *args)` command where we can optionally pass in specific initial conditions. This command will save these initial conditions inside the data folder and returns the foldername which will contain the data. 

We then run the `run_and_save` command which takes the initial conditions and `the_params` and executes the simulation. Under the hood, we check the boundary_conditions. If it is `periodic` (default) and (one-dimensional) we run `hkd1d()` to progress the simulation one timestep. Alternatively, we execute `hk1dreflective()` if it the boundary_conditions are reflective. `hegselmann_krause_update()` is used for 2d code. This code is stored in `_lib/model_helpers.jl`.

Once the simulation is complete, the output is saved to `data/$foldername/results.jld2`. This is then reloaded to manipulated to create the figures. The consensus statistics and the order parameters are defined in `_lib/plotting_helpers.jl`.

`phasediagrams.bash` executes a self-contained version of the code that has been optimised for efficiency due to the large cost of thousands of runs.

Random seeds are set at the beginning of each file and not for each simulation so the order of the code is relevant for reproducibility.