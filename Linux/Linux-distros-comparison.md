@def title="Linux distributions: identifying the ideal use cases."
@def date=2026-01-01
@def tag=["Linux"]

My name is Brenton Horne and I have been using Linux on and off since 2012, including several years in which I used various distributions as my daily driver. These distributions include, among others: Arch Linux, Debian, Fedora, Funtoo Linux, Gentoo Linux, Linux Mint, Mageia, Manjaro Linux, NixOS, OpenMandriva Lx, openSUSE, Sabayon and Ubuntu. Consequently, I would classify myself as an experienced user, and I wanted to give my opinion about the ideal use case of several Linux distributions, especially independent and innovative distributions. 

In the infoboxes I include in each distribution's section, I typically omit developmental releases when it comes to the release model and modernity sections. I do typically consider developmental releases when it comes to the initial release section. 

\tableofcontents 

~~~
<h1 style="clear: both;">Alpine Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="./Alpine_Linux_3.21.2.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">August 2005</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://alpinelinux.org/">alpinelinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">setup-alpine&mdash;text-based.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://wiki.alpinelinux.org/wiki/Alpine_Package_Keeper">Alpine Package Keeper</a> (APK; binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://wiki.alpinelinux.org/wiki/APKBUILD_Reference">APKBUILD</a>&mdash;shell script.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GNU Compiler Collection</a> (GCC)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/OpenRC">OpenRC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Musl">musl</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/BusyBox">BusyBox</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Almquist_shell">sh</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Detailed</td>
    </tr>
</table>
~~~

**Alpine Linux** is a security-focused distribution primarily intended for servers, routers, VPNs, and alike. A base Alpine Linux install can be as small as 144 MB in size and does not include Bash, sudo, Vim or nano.

~~~
<h1 style="clear: both;">Arch Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Arch_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">11 March 2002</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="http://www.archlinux.org/">www.archlinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer<exj/td>
        <td style="padding: 5px;"><a href="https://wiki.archlinux.org/title/Archinstall">archinstall</a>&mdash;textual installer.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://wiki.archlinux.org/title/Pacman">pacman</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://wiki.archlinux.org/title/PKGBUILD">PKGBUILD</a>&mdash;shell script.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a> (install), <a href="https://en.wikipedia.org/wiki/Z_shell">Zsh</a> (live).</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Vast, if the <a href="https://wiki.archlinux.org/title/Arch_User_Repository">Arch User Repository</a> (AUR) is included.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Comprehensive</td>
    </tr>
</table>
~~~
No conversation about Linux distributions geared towards advanced users would be complete without **Arch Linux**. It follows a Keep it Simple, Stupid (KISS) design philosophy. I may be biased in its favour as it is my go-to Linux distribution. A base install comes without a graphical user interface and has a pretty minimal array of packages, although the total size of a base install is about 1.7GB. It also has perhaps the most comprehensive documentation and vast repositories of any distribution.

It is ideal for users that:

* Are comfortable with the command line. Those not comfortable with the command line may favour EndeavourOS or Manjaro Linux.
* Want the freedom to customize their system, but without the desire to compile most components of their system from source. 
* Want the very latest software. On the flip side of this, they should also know how to recover from an update breaking their system.
* Prefer a rolling release model.
* Prefer a fast package manager. pacman is one of the fastest I have ever encountered. 
* May want obscure pieces of software. Packaging on Arch is easy for people familiar with shell script &mdash; the language of the Linux command line &mdash; and with its vast repositories many users do not even need to resort to packaging software for themselves. 
* Do not mind using standard system software like systemd. Users that dislike systemd may prefer Artix Linux. 

<h1 style="clear: both;">Chimera Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Chimera_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">2021</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://chimera-linux.org/">chimera-linux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Manual&mdash;bootstrapping.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://wiki.alpinelinux.org/wiki/Alpine_Package_Keeper">APK</a> (binary), <a href='https://github.com/chimera-linux/cports'>cports</a> (source)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://github.com/chimera-linux/cports">template.py</a>&mdash;Python script.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Clang">Clang</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://github.com/davmac314/dinit">Dinit</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Musl">musl</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/FreeBSD">FreeBSD</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Almquist_shell">sh</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
</table>

**Chimera Linux** (not to be confused with ChimeraOS) is a truly unique Linux distribution and uses a very unusual combination of system software components. One interesting characteristic of the distribution that I did not mention in the infobox to the right is that Chimera Linux does not come with `sudo` pre-installed and it does not seem to be in Chimera's repositories ([source](https://pkgs.chimera-linux.org/packages?name=sudo&origin=)). Given the distribution's security focus, as evidenced by its use of musl, I would imagine this omission is a deliberate security feature. 

The ideal use case of Chimera Linux would be on security-critical systems, with users that favour FreeBSD's command line, do not need vast repositories and prefer rolling release models. Especially those that prefer to write packages using Python scripts, prefer fast package managers and dislike systemd.

<h1 style="clear: both;">CRUX</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="CRUX_3.7.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">December 2002</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://crux.nu/">crux.nu</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Manual, with setup script.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">Ports with prt-get (source).</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://crux.nu/Main/Handbook3-7-Package">Pkgfile</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
</table>

**CRUX** aims to keep it simple as it uses tar.gz-based packages, BSD-style init scripts, and has fairly small repositories. It otherwise uses standard Linux system software. CRUX follows a fixed release model with new releases every year or two. It uses source-based package management and is best suited to advanced users that appreciate its idea of simplicity. A base install of CRUX 3.7, with GRUB installed to serve as the bootloader, uses about 2.6GB disk space. 

<h1 style="clear: both;">Debian</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Debian_12.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">August 1993</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://www.debian.org/">www.debian.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable<sup><a href="#footnote-5">5</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Debian-Installer">Debian-Installer</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">Advanced Packaging Tool</a> (APT; binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://wiki.debian.org/Packaging/Intro">Rules (Makefile), control, copyright and changelog files</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal or compete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Detailed</td>
    </tr>
</table>

**Debian** is the second-eldest Linux distribution still under active development. It has new stable releases every two years, roughly. The rules files used to build Debian packages are Makefiles with custom build commands. 

It has three and sometimes four branches. In ascending order of modernity, they are (fourth branch in brackets): (old stable), stable, testing and unstable. Old stable corresponds to the previous stable release of the distribution. The stable branch corresponds to the current stable release of the distribution; each stable release comes with three years of support. In the lead up to a new stable release, the testing branch is forked and frozen and the packages undergo further testing and potentially patching until they are ready to be incorporated in the next stable release. Unstable is where Debian's very latest packages start out, until after sufficient testing they make their way into testing. Testing and unstable branches follow a rolling release model and are cutting edge and bleeding edge, respectively. 

As users can choose a minimal install from its installer, and there are three main branches users can choose from, Debian can be a good choice for users that want to customize their system as much as one can without installing packages from source. Especially those that do not mind using systemd, like having very large repositories and do not mind having to use Makefiles to build packages, should one need to. Users needing more a beginner-friendly distribution should ideally use elementary OS, Linux Mint, MX Linux, Ubuntu and Zorin OS. 

<h1 style="clear: both;">Exherbo</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Exherbo.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">17 January 2006?<sup><a href="#footnote-6">6</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://exherbolinux.org/">exherbolinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Manual&mdash;bootstrapping and compiling.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://paludis.exherbolinux.org/">Paludis</a> (source)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">ebuild&mdash;shell script with custom commands.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a><sup><a href="#footnote-7">7</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**Exherbo** is a source-based distribution that originally was forked from Gentoo Linux. Like Gentoo, it uses ebuilds as its packaging files. Its package manager, Paludis, is written in C++ unlike Gentoo's Portage, which is written in Python. Paludis is specifically meant to be a better alternative to Portage. Given Exherbo has smaller repositories and less comprehensive documentation, but is practically the same as Gentoo except without Gentoo's binary repositories, I would be inclined to think that Exherbo is best suited to Gentoo fans that are disgruntled with Portage. 

<h1 style="clear: both;">Fedora</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Fedora_41.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">4 November 2003</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://fedoraproject.org/">fedoraproject.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Cutting edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Anaconda_(installer)">Anaconda</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/DNF_(software)">Dandified YUM</a> (DNF; binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**Fedora** is a distribution that can be argued to be beginner friendly, although I am inclined to not put it in that category because it does not have out-of-the-box support for proprietary drivers, including WiFi drivers. Fedora is one of the most up-to-date fixed release distributions I am aware of, although each release usually keeps to the same release (except for bug and security fix releases) of desktop environment software and with six months between releases, this makes it not truly bleeding edge. Fedora releasers come with 13 months of support, so users only need to upgrade to every second release, should they choose. Fedora also has an immutable root file system edition called Silverblue. Fedora is best suited to users that favour a fixed release model, like cutting edge software, need large repositories and prefer to package with spec files, when this is necessary. 

<h1 style="clear: both;">Gentoo Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Gentoo_Linux.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">31 March 2002</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://www.gentoo.org/">www.gentoo.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Cutting edge<sup><a href="#footnote-8">8</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Manual&mdash;bootstrapping and compiling.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Portage_(software)">Portage</a> (source)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">ebuild&mdash;shell script.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/OpenRC">OpenRC</a>/<a href="https://en.wikipedia.org/wiki/Systemd">systemd</a><sup><a href="#footnote-7">7</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Detailed</td>
    </tr>
</table>

**Gentoo Linux** is, in many ways, the prototypical source-based Linux distribution. It uses its own package manager called Portage which is meant to be a true ports system in the spirit of BSD ports. Although, in recent years it has become feasible to install most package as pre-compiled binaries. Interestingly, there have even been projects to port Portage to other operating systems like the BSD derivatives FreeBSD and NetBSD. It is ideal for Linux users that want complete freedom to customize their system all the way down to the configure/compile options used to build each package. Users can even fork packages and apply custom patches to them, should they choose.  

<h1 style="clear: both;">Guix System</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image src="Guix_System.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">29 March 2016<sup><a href="#footnote-9">9</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://guix.gnu.org/">guix.gnu.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Text-based installer.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Guix">GNU Guix</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">GNU Guile scripts.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Guix#GNU_Shepherd_init_system">GNU Shepherd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Detailed</td>
    </tr>
</table>

**Guix System** (pronounced "Geeks") is a reproducible, entirely free (as in freedom) Linux distribution that uses its own package manager called GNU Guix which installs each package to its own unique directory within `/gnu/store`. While Guix System uses GNU Guix as its package manager, GNU Guix is technically distribution-agnostic. Guix System is configured using files written in GNU Guile, such as `/etc/config.scm`. GNU Guile is also used to write packaging files for GNU Guix. Unlike NixOS, another reproducible Linux distribution, it does not seem to keep old configurations in its bootloader menu. It seems suitable for users that want a system entirely configurable using a single file written in GNU Guile and favour an entirely free operating system, even though this often comes with hardware compatibility issues. 

<h1 style="clear: both;">Linux From Scratch</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Linux_From_Scratch.png"/><caption style="text-align:left">The Linux From Scratch logo.</caption></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">3 December 1999</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="http://www.linuxfromscratch.org/">www.linuxfromscratch.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Manual compilation of each component.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">None, software manually compiled from source.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">None</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a>/<a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Detailed</td>
    </tr>
</table>

**Linux From Scratch** (LFS) is a source-based distribution wherein each software package is manually compiled and installed from source.

<h1 style="clear: both;">Mageia</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Mageia_9.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">1 July 2011</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://www.mageia.org/en">www.mageia.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">DrakX&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">DNF (current) and urpmi (legacy)&mdash;both binary.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Spec file.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**Mageia** is a Linux distribution that started out in 2011 as a fork of Mandriva Linux created by some former employees of the company that had developed Mandriva. Originally, it used the same tools as Mandriva like the package manager urpmi, but it has modernized in some ways and now uses DNF as its package manager. My experience with it is that it is rock solid stable, but many packages that I use are missing from their repositories. Consequently, I would recommend Mageia to users that want a rock solid stable system and do not have obscure software needs. Especially if they used Mandriva Linux and were fond of it.  

<h1 style="clear: both;">NixOS</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="NixOS_24.11.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">3 June 2003</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://nixos.org/">nixos.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">Nix (binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Nix file</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Comprehensive</td>
    </tr>
</table>

**NixOS** is a reproducible Linux distribution that is configured using files written in a special-made purely functional programming language called Nix. It also uses a package manager called Nix; Nix installs each package to its own unique directory within `/nix/store`. 

<h1 style="clear: both;">NuTyX</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="NuTyX_24.10.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">14 September 2009<sup><a href="#footnote-10">10</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://nutyx.org/en/">nutyx.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">NuTyX Installer&mdash;text-based.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">Cards (binary and source)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Pkgfile&mdash;shell script</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
</table>

**NuTyX** is based on LFS but has its own hybrid binary and source package manager called cards. 

<h1 style="clear: both;">openmamba GNU/Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="openmamba_GNU_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">Before 4 August 2009<sup><a href="#footnote-11">11</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://openmamba.org/">openmamba.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Cutting edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Mamba Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">DNF (binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
</table>

**openmamba GNU/Linux** is a Linux distribution that uses the DNF package manager. 

<h1 style="clear: both;">OpenMandriva Lx</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="OpenMandriva_Lx_Rolling.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">18 June 2013</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://www.openmandriva.org/">www.openmandriva.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed and rolling</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable (fixed), bleeding edge (rolling)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">urpmi (legacy) and DNF (current)&mdash;both binary.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Clang">Clang</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**OpenMandriva Lx** is a continuation of Mandriva Linux developed by a community project. Like Mageia, it uses the DNF package manager. Unlike Mageia, it comes in two editions&mdash;a fixed release and rolling release edition. It first started to offer a rolling release edition in 2023.<sup><a href="#footnote-12">12</a></sup> 

<h1 style="clear: both;">openSUSE</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="openSUSE_Tumbleweed.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">7 December 2006</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://www.opensuse.org/">www.opensuse.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed (Leap) and rolling (Tumbleweed)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable (Leap), bleeding edge (Tumbleweed)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">YaST&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">ZYpp (binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**openSUSE** is a continuation of the SUSE Linux distribution developed by a group of German computer science students and first released in March 1994. Like OpenMandriva Lx, it comes into separate editions &mdash; one, Leap, that features a fixed release model and another, Tumbleweed, that features a rolling release model. openSUSE adopted its two editions in 2014,<sup><a href="#footnote-13">13</a></sup> whereas OpenMandriva Lx adopted this two edition approach around 2023.<sup><a href="#footnote-12">12</a></sup>

<h1 style="clear: both;">PCLinuxOS</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="PCLinuxOS.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">October 2003</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://pclinuxos.com/">pclinuxos.com</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">PCLinuxOS Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/APT-RPM">APT-RPM</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**PCLinuxOS** is a beginner-friendly Linux distribution that was originally forked from Mandrake Linux 9.2 in 2003. It is rather conservative in some ways, for instance it still uses SysV as its init system, APT-RPM as its command-line package manager and Synaptic as its graphical package manager. APT-RPM had its last release in 2008, Synaptic has an outdated look although it is still maintained and SysV has been superseded on most distributions (not all, of course) by systemd, which was first released in 2010. Despite using a rolling release model, it also usually uses pretty old software. 

PCLinuxOS is perhaps best suited to beginners that do not need exotic software, like a no thrills and 2000s-esque desktop experience and favour a rolling release model. If somehow, despite being beginners, they have an opinion on init systems and dislike systemd, they may also like PCLinuxOS. 

<h1 style="clear: both;">Rhino Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Rhino_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">8 August 2023</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://rhinolinux.org/">rhinolinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://github.com/volitank/nala">Nala</a> (binary), <a href="https://github.com/pacstall/pacstall">Pacstall</a> (source) and distro-agnostic package managers (binary).<sup><a href="#footnote-14">14</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Per Debian, plus Pacscript (shell script).</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Vast</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**Rhino Linux** is a Ubuntu derivative that gets its core packages from the Ubuntu developmental branch. Rhino Linux is the only Ubuntu-based distribution that I am aware of that follows a rolling release model. It is specifically designed with developers in mind and comes with VSCodium pre-installed. It also comes with pacstall &mdash; a package manager that provides access to a repository designed to be the Ubuntu counterpart to the Arch User Repository &mdash; pre-installed. It also has a setup wizard that offers users four different distribution-agnostic package formats that the wizard can add support for onto their system &mdash; namely Snap, Flatpak, Nix and AppImages. Its default desktop environment is a customized Xfce desktop featuring a dock on the left of the screen that the distribution calls "Unicorn". 

<h1 style="clear: both;">Slackware Linux</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Slackware_Linux_15.0.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">17 July 1993</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="http://www.slackware.com/">www.slackware.com</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Slackware Installer&mdash;Textual.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">pkgtools</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://docs.slackware.com/slackware:slackbuild_scripts">SlackBuilds</a>&mdash;shell script</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

**Slackware Linux** is the eldest Linux distribution still in active development. A some unique characteristic of it is that its official repositories are fairly small and mostly just contain the packages that one can install from the live medium and updates thereof. There are unofficial repositories, but even they are not very large. This is largely because on Slackware it is expected that most non-core packages will be manually compiled from source using SlackBuild scripts. Another important characteristic of Slackware is that its developers are fairly conservative in that they are reluctant to adopt divisive pieces of technology like systemd, and often ship pretty old and well-tested versions of the software included in the system. I say they, but Slackware technically has a Benevolent Dictator for Life named Patrick Volkerding, who was its original creator back in 1993.

I have tried Slackware many times and I have found its approach to package management frustrating. That being said, it is rock solid and stable and if you are nostalgic for how Linux distros were like in the 1990s and like its approach to package management, it may be a suitable distribution for you. Especially if you are dislike the inclusion of systemd in most other modern Linux distributions.

<h1 style="clear: both;">Solus</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Solus.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">27 December 2016</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://getsol.us/">getsol.us</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;Textual.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Pardus_(operating_system)#PiSi_package_management">eopkg</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://help.getsol.us/docs/packaging/package.yml/">package.yml</a>&mdash;YAML file.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

Like PCLinuxOS, **Solus** is a Linux distribution aimed towards beginners, despite following a rolling release model. As its package manager, it uses eopkg&mdash;which is based on <a href="https://en.wikipedia.org/wiki/Pardus_(operating_system)">Pardus</a>' abandoned package manager of PiSi. It is the original distribution that the Budgie desktop was developed for, and is noted for its relatively good default aesthetics. 

Its ideal use case is probably a beginner that does not need exotic software, appreciates a beautiful and simple desktop like Budgie and does not want to have to upgrade their system between releases of a fixed release distribution like Linux Mint for fear of system breakage.

<h1 style="clear: both;">Vanilla OS</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Vanilla_OS_2.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">29 December 2022<sup><a href="#footnote-15">15</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://vanillaos.org/">vanillaos.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Vanilla Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;">ostree (read-only root), distro-agnostic package managers<sup><a href="#footnote-16">16</a></sup> and <a href="https://github.com/Vanilla-OS/apx">Apx</a> (applications).</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;">Any, due to Apx.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a><td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Vast</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Medium</td>
    </tr>
</table>

I would like to start this section by mentioning that I have not been able to install **Vanilla OS**, as the installer repeatedly fails for me. That being said, Vanilla OS uses a Debian (unstable) base (although, it previously used a Ubuntu base) and an immutable root file system. Actually, it has two root file systems. One is booted by the user, the other is the one to which updates are applied. This is so that users can boot this updated system at their next reboot but have their other root file system as a backup should the update break their system. Vanilla OS uses Apx (pronounced "apex") which uses containerized Linux distributions to provide access to software packaged for that distribution. 

<h1 style="clear: both;">Void</h1>
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img src="Void.png"/></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Initial release</td>
        <td style="padding: 5px;">2008</td>
    </tr>
    <tr>
        <td style="padding: 5px; width: 150px;">Website</td>
        <td style="padding: 5px;"><a href="https://voidlinux.org/">voidlinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Release model</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 5px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Installer</td>
        <td style="padding: 5px;">Void Installer&mdash;textual.</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Package manager (type)</td>
        <td style="padding: 5px;"><a href="https://github.com/void-linux/xbps">X Binary Package System</a> (XBPS; binary)</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Packaging file(s)</td>
        <td style="padding: 5px;"><a href="https://xbps-src-tutorials.github.io/">template</a>&mdash;shell script</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Compiler</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Init system</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Runit">runit</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">C standard library</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a>/<a href="https://en.wikipedia.org/wiki/Musl">musl</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Userland</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Shell</td>
        <td style="padding: 5px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 5px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 5px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 5px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 5px;">Minimal</td>
    </tr>
</table>

**Void** is a Linux distribution that is similar to Arch Linux in that it follows a rolling release model and uses a fast, lightweight package manager written in C that uses shell script packaging files for building its packages. A base install without a GUI is about 3GB in size.

<h1 style="clear: both;">Footnotes</h1>
<ol>
    <li id="footnote-1">The possible categories are, in ascending order of modernity, old stable, stable, cutting edge and bleeding edge.</li>
    <li id="footnote-2">It is difficult to rate this in a completely non-subjective way. As while you can list the number of packages in their repositories, some distributions package split a single piece of software into multiple packages and hence raw numbers are not as fair a measure as they seem. Plus some distributions have multiple variants on more or less the same package in their repositories. To simplify things, I will categorize repository size as: vast, very large, large, medium, medium-small, small and tiny.</li>
    <li id="footnote-3">A distribution will have "minimal" in this category if the base installation comes without a graphical user interface (GUI). If a GUI is a required part of a base install, I will say "complete".</li>
    <li id="footnote-4">The categories, in ascending order of documentation sufficiency, are: minimal, medium, detailed and comprehensive. Minimal documentation would typically just cover installation, package management and basic configuration. Medium would cover more configuration options than just basic. Detailed would cover some additional topics. Comprehensive would cover almost every conceivable topic for the distribution.</li>
    <li id="footnote-5">This is for the stable branch. Testing and unstable &mdash; which are the distribution's developmental branches &mdash; are cutting edge and bleeding edge, respectively.</li>
    <li id="footnote-6">I see that Paludis, Exherbo's package manager, has its first <a href="https://paludis.exherbolinux.org/changelog.html">changelog</a> entry on 17 January 2006. This is why I assume this was when Exherbo was established.</li>
    <li id="footnote-7">This is the default, but there are instructions for how to install other init systems.</li>
    <li id="footnote-8">This is the default. It is possible to install more bleeding edge software by adjusting Portage's keyword settings.</li>
    <li id="footnote-9">The earliest news about Guix System that I found was that the <a href="https://web.archive.org/web/20150203220723/http://www.fsf.org/news/fsf-adds-guix-system-distribution-to-list-of-endorsed-distributions">Free Software Foundation added it to their list of free Linux distributions on 3 February 2015</a>. The earliest blog post on the GNU Guix website that mentions Guix System was from <a href="https://guix.gnu.org/en/blog/2015/gnu-guix-talk-at-opentechsummit-berlin-may-14th/">12 May 2015</a>. The <a href="https://guix.gnu.org/en/blog/2016/gnu-guix--guixsd-0100-released/">first announcement</a> of a release that I found on their website was from 29 March 2016.</li>
    <li id="footnote-10">Earliest release in DistroWatch's database is the 2009 release which was <a href="https://distrowatch.com/5667">reported on 14 September 2009</a>. NuTyX's website has a copyright notice that begins in 2007, so its first release could be as long ago as 2007.</li>
    <li id="footnote-11">Earliest release in DistroWatch's database is the 1.1 release which was <a href="https://distrowatch.com/?newsid=05607">released on 4 August 2009</a>. This was meant to be an update to the earlier release of openmamba GNu/Linux 1.0.</li>
    <li id="footnote-12">This I say based on <a href="https://distrowatch.com/?newsid=11735">this news release from DistroWatch</a>.</li>
    <li id="footnote-13">Source: <a href="https://en.opensuse.org/Portal%3ATumbleweed?">Portal:Tumbleweed at openSUSE Wiki</a>.</li> <!-- openSUSE-->
    <li id="footnote-14">I mention these, even though I omit them in most distro's infoboxes, because Rhino Linux has options to enable cross-distro package managers/formats in its setup wizard. Specifically, it allows users to enable Flatpak, Nix, Snap or AppImages.</li> 
    <li id="footnote-15"><a href="https://vanillaos.org/blog/article/2022-12-29/vanilla-os-2210-kinetic-is-out">Vanilla OS 22.10</a> was this release and it was the first release mentioned in the <a href="https://vanillaos.org/blog">Vanilla OS</a> blog.</li>
    <li id="footnote-16">I say this because <a ref="https://docs.vanillaos.org/handbook/en/install-flatpaks">Flatpak installation instructions</a> and <a href="https://docs.vanillaos.org/handbook/en/install-homebrew">Homebrew installation instructions</a> are given in Vanilla OS's documentation.</li>
</ol>
