# ICA-PCA
This code seperates 3 mixed signals to their original sources using ICA.
  1. 3 mixed sound files without noise (mix1-3) are seperated easily with ICA
  2. 3 mixed sound files with an added noise (noisy_mix1-3) are attempted to be seperated with ICA. the result is poor
  3. the sound files with noise are attempted to be seperated, this time PCA is applied on 7 recordings to reduce the dimensions.
  then ICA is applied to seperate the mixed recordings to the same 3 original sources. the results are much better.
