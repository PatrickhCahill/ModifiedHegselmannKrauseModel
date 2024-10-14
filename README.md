# ModifiedHegselmannKrauseModel
This repository is designed to reproduce the simulated figures from the *A modified Hegselmann–Krause model for interacting voters and political parties* (2024) by Patrick Honnery Cahill and Georg Gottwald.

To reproduce the figures [juliaLANG needs to be installed](https://julialang.org/downloads/) - bash will also be required. There are four separate files that will need to be executed. First, open a terminal and execute `mkdir data` (which is automatically gitignored) and then execute the following either in the terminal or from the REPL:
* `requirements.kl` This will install the necessary packages in order to execute the code.
* `voter_behaviours_explanation.jl` This will create each subfigure in Figure 2 from the paper.
* `figure3.jl` will create figure 3. **Manual editing to point to the correct foldername is required.**
* `4to8pipeline.jl` This populate the figures from 4 to 8.
* `phasediagrams.bash` will create execute the heatmap related julia files in parallel in order to create figures 10 and 11. **Note requires bash and is computationally expensive.**


Note: Figures are provided with a seed to be reproducible. But they may differ between versions of julia and package versions. The figures in text were produced using `julia version 1.11.0-beta1`.