Summary: Installer for Beagle
License:                GPL
Name:                   beagle
Version:                2.1
Release:                0
Group:                  Development/Tools
Source: beagle-2.1.tar.gz

%description
Beagle provides performance boost for packages like BEAST, MrBayes, and GARLI.

%prep

%setup -q

%build
./autogen.sh
./configure --prefix=/usr/local --disable-rpath
make

%install
rm -rf $RPM_BUILD_ROOT
make prefix=$RPM_BUILD_ROOT/usr/local install

%clean
rm -rf $RPM_BUILD_ROOT

%files
/usr/local/include/libhmsbeagle*
/usr/local/lib/libhmsbeagle*
/usr/local/lib/pkgconfig/hmsbeagle*

%defattr(-,root,root)
