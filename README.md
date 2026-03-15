Moiree Spectrum Analizer - a research and entertainment tool to look at LFM-CW radar signal forms in Fourier space depending on changes in various parameters of the signal model.

Kinda, let's just wait for a paper on my website. It should be much cooler there with illustrated letters.

## Build: ##
```
cmake -S . -B build
cd build
make
make install
```

The program and (probably) all necessary dependencies will reside in the "out" directory in the source root. If the program crashes at startup, it's probably due to lack of "data" directory inside the "out" one. Copy "data" directory from the source root to "out" directory -- and the program should start properly.