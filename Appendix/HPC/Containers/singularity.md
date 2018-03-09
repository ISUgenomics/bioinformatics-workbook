#What is Singularity?


##How does Singularity overcome common issues with containers?

##How do I use Singularity?


## Initial setup

The first time you use singularity it will by default put a .singularity folder in your home directory which commonly has limited storage space.  Therefore it is important that you move that folder to a different location and then create a softlink from your home directory to the new location.


## Where to find containers.

* https://singularity-hub.org/
* https://hub.docker.com


## How to create a singularity image from a docker image


```
singularity pull docker://sjackman/maker
singularity exec --bind $PWD ./maker.img maker --help
```
