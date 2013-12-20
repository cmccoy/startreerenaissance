
place in /root/rpmbuild

  rsync -rv ./ root@remote:rpmbuild

install deps

  yum install -y rpm-build autoconf automake tree libtool
