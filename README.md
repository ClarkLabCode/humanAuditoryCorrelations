# humanAuditoryCorrelations
Code for generating stimuli and analyzing data related to human auditory correlation sensitivity.

## Each folder corresponds to a psychophysical task:

**ternaryExperiment**: contains the psychophysical script (`ternary.m`), tone generating function (`ternaryFunction.m`), and analysis code (`analyzeTernaryExperiment.m`) to generate **Figure 1c**. 

**coherenceExperiment**: contains the psychophysical script (`coherence.m`), tone generating function (`coherenceFunction.m`), and analysis code (`analyzeCoherenceExperiment.m`) to generate **Figure 1d**. 

**binauralExperiment**: contains the psychophysical script (`binaural.m`), tone generating function (`binauralFunction.m`), and analysis code (`analyzeBinauralExperiment.m`) to generate **Figure 1f**. 

**pipDeltaTimeExperiment**: contains the psychophysical script (`pipDeltaTime.m`), tone generating function (`pipDeltaTimeFunction.m`), and analysis code (`analyzePipDeltaTimeExperiment.m`) to generate **Figure 2c**. 

**pipDeltaTimeSuppExperiment**: contains the psychophysical script (`pipDeltaTimeSupp.m`), tone generating function (`pipDeltaTimeSuppFunction.m`), and analysis code (`analyzePipDeltaTimeSuppExperiment.m`) to generate **Figure S2a, b**. 

**pipDeltaNoteExperiment**: contains the psychophysical script (`pipDeltaNote.m`), tone generating function (`pipDeltaNoteFunction.m`), and analysis code (`analyzePipDeltaNoteExperiment.m`) to generate **Figure 2d**. 

**correlatedPipsExperiment**: contains the psychophysical script (`correlatedPips.m`), tone generating function (`correlatedPipsFunction.m`), and analysis code (`analyzeCorrelatedPipsExperiment.m`) to generate **Figure 2b, c**. 

**glidersExperiment**: contains the psychophysical script (`gliders.m`), tone generating function (`glidersFunction.m`), and analysis code (`analyzeGlidersExperiment.m`) to generate **Figure S3c, d**.

**soundAnalysisCode**: [TO-DO]

**simulations**: contains scripts `simulateTrackingHeuristics.m`, `simulateMotionEnergyUnitResponsesSweep.m`, and `simulateMotionEnergyRespToAllStim.m` to generate **Figures S2c, d and S3a; 4a, b, c; S4a, b, c**, respectively. `makeRBcolormap.m` and `niceAxesLarge.m` are formatting helper functions that need to be in the working folder in order to run the scripts.

**generalToneGeneration**: contains the tone generating function `generalToneGenerationFunction.m` for creating and listening to the tone types that were used in these experiments.

## Specifications:
All psychophysical task code was run on a MacBook Pro with an Intel chip using Matlab 2021b and Psychtoolbox.

## Directions:

1) Ensure that the above specifications are met.

2) To participate in a task, download the psychophysical script and corresponding tone generating function. Make sure that `iso226.m` and the corresponding data folder are also in the working folder. Run the psychophysical script, and enter a subject name of choice. Follow the task instructions until it concludes.
   
3) To analyze psychophysical data from the original experiments, use Dryad [INSERT LINK HERE; TO-DO] to download the corresponding data folder for each experiment, which contains the `.mat` files for all the participants for a particular experiment. Make sure that the analysis script is in the same folder as the data before running the script.
