Codes to extract paddle-to-paddle corrections:
1) parameter_extraction.cpp -> do the extraction by using the bandHits that have already been corrected for time-walk and L-R.

2) parameter_extraction_rawBanks.cpp -> do the extraction on the raw banks and apply corrections by hand.

*************************************************************
                  ___                  ___
\    /\    / /\   |  \ |\  | | |\  |  /
 \  /  \  / /--\  |-\/ | \ | | | \ | |  --\
  \/    \/ /    \ |  \ |  \| | |  \|  \___/

Run this code on laser data

If you are going to use 'parameter_extraction_rawBanks.cpp',
make sure you have the most updated CCDB-style files needed
in 'input_params'

*************************************************************

We align all the bars in each layer to a reference bar in that layer. Then, as a second step, we align each layer to a reference layer.
To be more specific, we align all the bars from layer i to bar corresponding to (sector: 2, layer: i, component: 1).
Then, we align these reference bars (sector: 2, layer: i, component: 1) for each layer to bar (sector: 2, layer: 5, component: 1).
