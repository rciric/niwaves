# niwaves

Translating lag (and possibly complexity eventually)

## What is `niwaves`?

Many tools used to interrogate complex systems in neuroimaging are written in Matlab. `niwaves` initially aimed either to translate these tools or to index existing translations in order to open them to the Python ecosystem. The scope of this project is rather limited, so we hope that our efforts are cannibalised into a more mature toolbox at some point in the future.

As of now, this library contains scripts for lagged correlation analysis following Mitra et al. The code needs to be tested more extensively; for now, if you plan on using this, we strongly suggest that you cross-check results against the existing Matlab version before proceeding.

On our to-do list (maybe):

1. Lag and lag threads (sketched as of December 2018)

2. [Network community toolbox](http://commdetect.weebly.com): translation of functions by Dr. Shi Gu available [here](https://github.com/nangongwubu/Python-Version-for-Network-Community-Architecture-Toobox)

3. [Network control](https://www.danisbassett.com/research-projects.html), translation in progress elsewhere

4. [Algebraic topology](https://www.aesizemore.com/network-toolboxes.html)

Central to implementing many of these tools will be operationalising a substitute for the excellent [`genlouvain`](http://netwiki.amath.unc.edu/GenLouvain/GenLouvain). There's a [genlouvain function included](https://github.com/nangongwubu/Python-Version-for-Network-Community-Architecture-Toobox/blob/master/build/lib/ncat.py#L916) in Dr. Shi Gu's Network Community Architecture Toolbox.
