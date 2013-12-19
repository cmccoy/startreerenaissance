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
./configure --prefix=%{_prefix}
make
 
%install
rm -rf $RPM_BUILD_ROOT
make prefix=$RPM_BUILD_ROOT%{_prefix} install
 
%clean
rm -rf $RPM_BUILD_ROOT
 
%files
%{_prefix}/include/libhmsbeagle*
%{_prefix}/lib/libhmsbeagle*
%{_prefix}/lib/pkgconfig/hmsbeagle*
 
%defattr(-,root,root)
