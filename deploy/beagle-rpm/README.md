# README

Files in this directory assist in building an RPM for [BEAGLE](https://code.google.com/p/beagle-lib/).

## Steps to build on EC2 amazon linux node:

### place in `/root/rpmbuild`

    rsync -rv ./ root@remote:rpmbuild

###  On remote systeM:

#### install deps

    yum install -y rpm-build autoconf automake libtool

#### build

    cd rpmbuild
    rpmbuild -bc SPECS/beagle.spec
