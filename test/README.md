markovtest
==========

This is for running a docker container that downloads and makes the sample code from the Semi-Markov library.

```
sudo docker run u-gcc48-boost55 /bin/bash -c "git clone https://github.com/adolgert/markovtest.git; cd markovtest; ./build.sh"
```

And to make the images

```
for i in u-*; do docker build -t $i - < $i; done
```
