# Josephson Junction optimization via enhanced Shortcuts to Adiabaticity
In this repo I will try to use [eSTA](https://arxiv.org/abs/1912.06057) to minimize the decoherence of the final state of the system in the context of the Josephson Junction.
The files are organized as follows:
```

.
├── fidelityplot.jl
├── Manifest.toml
├── Project.toml
├── README.md
├── src
│   ├── ConstantQuantities.jl
│   ├── Fidelities.jl
│   ├── ModulationNoise.jl
│   ├── PlotObj.jl
│   ├── Squeezing.jl
│   ├── TimeNoise.jl
│   └── WhiteNoise.jl
└── test
    └── whitenoise.jl
```
The `src` folder contains the main functions used in the calculations, while the `test` folder contains the scripts that are used to test the functions in the `src` folder.
This calculations rely on the [JosephsonSTA package](github.com/manueloid/JosephsonSTA) which is a Julia package that implements the eSTA algorithm.

I have three different types of noise for the moment
- Modulation noise: of the form Λδ(t) = (1+δ) Λ(t)
- Time shift noise: of the form Λδ(t) = Λ(t+δ)
- White noise: of the form Λδ(t) = Λ(ξ(t)), where ξ(t) is a Gaussian noise

For each of them, I will put the relevant scripts in the `src` folder.
I will name the scripts as follows:
- `ModulationNoise.jl`
- `TimeShiftNoise.jl`
- `WhiteNoise.jl`

I will use the camel case for the files that contain the functions that are used in the calculations, and the  lowercase for the files that contain the scripts that are used to calculate the fidelities and the robustness of the system.
