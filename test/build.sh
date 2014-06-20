#!/bin/bash
export PATH=/usr/bin:$PATH
mkdir /home/jenkins
cd /home/jenkins
git clone https://github.com/afidd/Semi-Markov.git
cd Semi-Markov
echo "autoreconf now"
autoreconf --force --install
echo "done autoreconf"
./configure
make
./sirmixed
./meta
./weiss
./bvd
