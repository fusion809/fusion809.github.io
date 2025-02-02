@def title="Linux distributions: identifying the ideal use cases."
@def date=2026-01-01
@def tag=["Linux"]

My name is Brenton Horne and I have been using Linux on and off since 2012, including several years in which I used various distributions as my daily driver. These distributions include, among others: Arch Linux, Debian, Fedora, Funtoo Linux, Gentoo Linux, Linux Mint, Mageia, Manjaro Linux, NixOS, OpenMandriva Lx, openSUSE, Sabayon and Ubuntu. Consequently, I would classify myself as an experienced user, and I wanted to give my opinion about the ideal use case of several Linux distributions, especially independent and innovative distributions. 

In the infoboxes I include in each distribution's section, I typically omit developmental releases when it comes to the release model and modernity sections. I do typically consider developmental releases when it comes to the initial release section, however. The images I show are largely hyfetch, neofetch or fastfetch output. For Linux From Scratch, I just used the official logo.

When I mention "exotic" or "obscure" software, I mean software that is fairly unpopular and used for niche purposes. For instance, the Marvin Suite of ChemAxon is a piece of software for sketching skeletal formulas, among other things, and I would class it as exotic or obscure as it is only used for a fairly niche purpose and not many people are aware of it.

~~~
<table style="float: left; border-collapse: collapse;">
<tr>
    <td style="font-size: 30px; padding: 10px;"><b>Table of contents</b></td>
</tr>
<tr>
    <td style="font-size: 20px; padding: 10px;">
    <ol>
        <li><a href="#preliminaries">Preliminaries</a></li>
        <ul>
        <li><a href="#what_is_linux">What is Linux?</li>
        <li><a href="#core_components_of_a_Linux_distribution">Core components of a Linux distribution</a></li>
        <li><a href="#linux_graphical_user_interface">Linux graphical user interface</a></li>
        <li><a href="#release_model">Release model</a></li>
        </ul>
        <li><a href="#alpine_linux">Alpine Linux</a></li>
        <li><a href="#arch_linux">Arch Linux</a></li>
        <li><a href="#chimera_linux">Chimera Linux</a></li>
        <li><a href="#crux">CRUX</a></li>
        <li><a href="#debian">Debian</a></li>
        <li><a href="#deepin">deepin</a></li>
        <li><a href="#elementary_os">elementary OS</a></li>
        <li><a href="#exherbo">Exherbo</a></li>
        <li><a href="#fedora">Fedora</a></li>
        <li><a href="#gentoo_linux">Gentoo Linux</a></li>
        <li><a href="#guix_system">Guix System</a></li>
        <li><a href="#linux_from_scratch">Linux From Scratch</a></li>
        <li><a href="#linux_mint">Linux Mint</a></li>
        <li><a href="#mageia">Mageia</a></li>
        <li><a href="#mx_linux">MX Linux</a></li>
        <li><a href="#nixos">NixOS</a></li>
        <li><a href="#nutyx">NuTyX</a></li>
        <li><a href="#openmamba_gnulinux">openmamba GNU/Linux</a></li>
        <li><a href="#openmandriva_lx">OpenMandriva Lx</a></li>
        <li><a href="#opensuse">openSUSE</a></li>
        <li><a href="#pclinuxos">PCLinuxOS</a></li>
        <li><a href="#rhino_linux">Rhino Linux</a></li>
        <li><a href="#slackware_linux">Slackware Linux</a></li>
        <li><a href="#solus">Solus</a></li>
        <li><a href="#ubuntu">Ubuntu</a></li>
        <li><a href="#vanilla_os">Vanilla OS</a></li>
        <li><a href="#void">Void</a></li>
        <li><a href="#footnotes">Footnotes</a></li>
    </ol>
    </td>
</tr>
</table>
~~~

# Preliminaries
## What is Linux?
A Linux distribution is defined by its use of the [Linux kernel](https://en.wikipedia.org/wiki/Linux_kernel). The kernel is the component of an operating system that manages all communication between hardware and software. [File system](https://en.wikipedia.org/wiki/File_system) (a system wherein files are stored) support is also within the kernel, although some tools for managing file systems (e.g. `mkfs` command) are provided by other core system software. Popular file systems with Linux kernel support include [Btrfs](https://en.wikipedia.org/wiki/Btrfs), [ext4](https://en.wikipedia.org/wiki/Ext4), [FAT32](https://en.wikipedia.org/wiki/File_Allocation_Table#FAT32) (popular for [EFI system partitions](https://en.wikipedia.org/wiki/EFI_system_partition)) and [XFS](https://en.wikipedia.org/wiki/Xfs). [ZFS](https://en.wikipedia.org/wiki/ZFS) also has support on Linux via a third-party kernel module. This is because ZFS' licence is incompatible with the kernel's [GNU General Public Licence](https://en.wikipedia.org/wiki/GNU_General_Public_License) (GPL) and hence cannot be included in the kernel's source code.

Linux distributions are almost always [Unix-like](https://en.wikipedia.org/wiki/Unix-like) too, although there are exceptions like [Android](https://en.wikipedia.org/wiki/Android_(operating_system)). This means, among other things, that most Linux distributions share similar commands to Unix systems and roughly follow the design philosophy of Unix with each command doing just one thing and doing it well. 

## Core components of a Linux distribution
What other software a Linux distribution uses varies markedly between distributions and even between different installs of the same distribution. Other core components of a Linux operating system include, but are not limited to:
* [Init system and service manager](https://en.wikipedia.org/wiki/Init). The init is the first process that the kernel starts after the system is booted. It starts, either directly or indirectly, all other processes. Service managers manage [daemons](https://en.wikipedia.org/wiki/Daemon_(computing)), which are background processes. These processes can be important and manage things like the internet connection (e.g. dhcpcd) and graphical login manager. The most popular Linux init systems, many of which also function as service managers, are in approximately descending order of popularity: [systemd](https://en.wikipedia.org/wiki/Systemd), [SysV init](https://en.wikipedia.org/wiki/Init#SysV-style), [OpenRC](https://en.wikipedia.org/wiki/OpenRC), [runit](https://en.wikipedia.org/wiki/Runit), [Dinit](https://davmac.org/projects/dinit/), [BusyBox-init](https://en.wikipedia.org/wiki/BusyBox), [GNU Shepherd](https://en.wikipedia.org/wiki/GNU_Guix#GNU_Shepherd_init_system) and [s6](https://skarnet.org/software/s6/). systemd is disliked by some Linux users, one reason they often cite for this is that it violates the Unix philosophy in that it does so many different things (including providing its own replacements for some daemons). 
* [C standard library](https://en.wikipedia.org/wiki/C_standard_library). This provides functions, macros and type definitions that programs written in C&mdash;which many core components of any Linux system are&mdash;can call. The vast majority of Linux distributions used [glibc](https://en.wikipedia.org/wiki/Glibc) (from the [GNU Project](https://en.wikipedia.org/wiki/GNU_Project)) to provide this. [Musl](https://en.wikipedia.org/wiki/Musl) is a less commonly used alternative that has a security focus. [uClibc](https://en.wikipedia.org/wiki/UClibc) is another possible alternative for Linux distributions.
* [Toolchain, including a compiler](https://en.wikipedia.org/wiki/Toolchain). This is basically what it used to build most components of the system, converting source code into binaries that can be executed. [GNU Toolchain](https://en.wikipedia.org/wiki/GNU_toolchain), including the [GNU Compiler Collection](https://en.wikipedia.org/wiki/GNU_Compiler_Collection) (GCC), is the most common toolchain on Linux. [LLVM](https://en.wikipedia.org/wiki/LLVM) and [Clang](https://en.wikipedia.org/wiki/Clang) are popular alternatives.
* [The Unix shell](https://en.wikipedia.org/wiki/Unix_shell). This serves as a crucial part of the systems command-line interface. [Bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) is the most common Unix shell on Linux and is developed as part of the GNU Project. [Zsh](https://en.wikipedia.org/wiki/Z_shell) is a popular alternative, although it is not a popular default for Linux distributions. [tcsh](https://en.wikipedia.org/wiki/Tcsh) and [sh](https://en.wikipedia.org/wiki/Almquist_shell) are less popular Unix shells for Linux. 
* [Unix commands, or userland utilities](https://en.wikipedia.org/wiki/List_of_POSIX_commands). These are used to perform common command-line tasks such as copying files, moving files, manipulating strings, producing checksums, determining the current directory, etc. Most Linux distributions use the [GNU coreutils](https://en.wikipedia.org/wiki/GNU_Core_Utilities) package to provide these commands. [BusyBox](https://en.wikipedia.org/wiki/BusyBox) also provides many of these commands. Other Unix and Unix-like systems such as [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD), [NetBSD](https://en.wikipedia.org/wiki/NetBSD), [OpenBSD](https://en.wikipedia.org/wiki/OpenBSD) and [OpenIndiana](https://en.wikipedia.org/wiki/OpenIndiana) provide their own versions of these commands. I know that at least FreeBSD's Unix commands have been ported to Linux.
* [util-linux](https://en.wikipedia.org/wiki/Util-linux). It is provided by the Linux Kernel Organization, like the kernel itself. util-linux provides some more commands on most Linux distributions such as `chsh` for changing the default Unix shell, `dmesg` for checking the status of kernel modules and `fsck` for checking and fixing file systems. 
* [Package manager](https://en.wikipedia.org/wiki/Package_manager). This provides a means of installing, removing and updating software on one's system. There are several different package managers, including, but not limited to:
    - [dpkg](https://en.wikipedia.org/wiki/Dpkg) and its frontend [Advanced Packaging Tool](https://en.wikipedia.org/wiki/APT_(software)) (APT).
    - [RPM](https://en.wikipedia.org/wiki/RPM_Package_Manager) and its frontends [APT-RPM](https://en.wikipedia.org/wiki/APT-RPM), [Dandified YUM](https://en.wikipedia.org/wiki/Dandified_YUM) (DNF), [urpmi](https://en.wikipedia.org/wiki/Urpmi), [YUM](https://en.wikipedia.org/wiki/Yum_(software)) and [ZYpp](https://en.wikipedia.org/wiki/ZYpp).
    - [pacman](https://en.wikipedia.org/wiki/Arch_Linux#Pacman).
    - [Portage](https://en.wikipedia.org/wiki/Portage_(software)).

As many core components of most Linux distributions are created as part of the GNU Project, many argue that most Linux distributions should be called GNU/Linux distributions to give appropriate credit to the GNU Project.

## Linux graphical user interface
Most desktop Linux installations have a [graphical user interface](https://en.wikipedia.org/wiki/Graphical_user_interface) (GUI) too. A sufficiently complete Linux GUI that implements the desktop metaphor is often called a [desktop environment](https://en.wikipedia.org/wiki/Desktop_environment). Components of Linux GUIs include:
* [Graphics libraries](https://en.wikipedia.org/wiki/Graphics_library), such as [Mesa](https://en.wikipedia.org/wiki/Mesa_(computer_graphics)). These provide optimized functions for other GUI components. 
* [Display server communication protocol](https://en.wikipedia.org/wiki/Windowing_system#Display_server_communications_protocols), such as [X.Org Server](https://en.wikipedia.org/wiki/X.Org_Server), [Wayland](https://en.wikipedia.org/wiki/Wayland_(protocol)) and [Mir](https://en.wikipedia.org/wiki/Mir_(software)). X.Org Server has been the dominant protocol on Linux and other Unix/Unix-like systems since the 1980s, but is in the process of being replaced by Wayland. Mir was originally developed with the intention of it replacing X.Org Server on the desktop, but it now seems to be largely developed for use in the Internet of Things. 
* [Widget toolkit](https://en.wikipedia.org/wiki/Widget_toolkit) to build GUIs for programs. Examples include [Enlightenment Foundation Libraries](https://en.wikipedia.org/wiki/Enlightenment_Foundation_Libraries) (EFL), [GTK](https://en.wikipedia.org/wiki/GTK), [Qt](https://en.wikipedia.org/wiki/Qt_(software)), [Tk](https://en.wikipedia.org/wiki/Tk_(software)), etc.
* [X window manager](https://en.wikipedia.org/wiki/Window_manager) (for X.Org Server)/[Wayland compositor](https://en.wikipedia.org/wiki/Wayland_compositor). This manages the windows for the GUI. Examples include [Mutter](https://en.wikipedia.org/wiki/Mutter_(software)), [KWin](https://en.wikipedia.org/wiki/KWin), [Enlightenment](https://en.wikipedia.org/wiki/Enlightenment_(window_manager)), [i3](https://en.wikipedia.org/wiki/I3_(window_manager)), [Sway](https://en.wikipedia.org/wiki/Sway_(window_manager)) and [Muffin](https://en.wikipedia.org/wiki/Muffin_(software)).
* [Graphical shell](https://en.wikipedia.org/wiki/Shell_(computing)#Graphical_shells), which provides menus, panels, docks, taskbars, system notifications, desktop icons, wallpaper management, application launchers, etc. Examples include the [GNOME Shell](https://en.wikipedia.org/wiki/GNOME_Shell), [Unity](https://en.wikipedia.org/wiki/Unity_(user_interface)) and that of [KDE Plasma](https://en.wikipedia.org/wiki/KDE_Plasma).
* [Graphical login manager](https://en.wikipedia.org/wiki/X_display_manager) (also known as a display manager), such as [GNOME Display Manager](https://en.wikipedia.org/wiki/GNOME_Display_Manager) (GDM), [LightDM](https://en.wikipedia.org/wiki/LightDM) and [Simple Desktop Display Manager](https://en.wikipedia.org/wiki/Simple_Desktop_Display_Manager) (SDDM).

Complete desktop environments include [GNOME](https://en.wikipedia.org/wiki/GNOME) (originally part of the GNU Project), KDE Plasma (versions [4](https://en.wikipedia.org/wiki/KDE_Plasma_4), [5](https://en.wikipedia.org/wiki/KDE_Plasma_5) and [6](https://en.wikipedia.org/wiki/KDE_Plasma_6)), [Deepin Desktop Environment](https://en.wikipedia.org/wiki/Deepin#Deepin_Desktop_Environment), [Cinnamon](https://en.wikipedia.org/wiki/Cinnamon_(desktop_environment)), [Budgie](https://en.wikipedia.org/wiki/Budgie_(desktop_environment)), [LXDE](https://en.wikipedia.org/wiki/LXDE), [MATE](https://en.wikipedia.org/wiki/MATE_(desktop_environment)), [LXQt](https://en.wikipedia.org/wiki/LXQt), [Trinity](https://en.wikipedia.org/wiki/Trinity_Desktop_Environment), [UKUI](https://en.wikipedia.org/wiki/UKUI), [Lumina](https://en.wikipedia.org/wiki/Lumina_(desktop_environment)) and [COSMIC](https://en.wikipedia.org/wiki/COSMIC_(desktop_environment)). 

## Release model
Another difference between Linux distributions is their release model. Release models come in three different categories: fixed, rolling and semi-rolling. 

Most Linux distributions follow a fixed release model with new releases of the distribution every so often. Typically, these releases differ in the versions of system software and GUI software included in the system. These pieces of system software do often receive security and bug fixes over the lifetime of a given release, but they usually do not undergo major updates during that time. Most distributions following a fixed release have a mechanism to upgrade from one release to the next. Although, these upgrades can cause system breakage. 

Some distributions follow a [rolling release model](https://en.wikipedia.org/wiki/Rolling_release) in which there are no fixed releases of the system. Rather updates to each system and GUI component are just provided whenever they are ready, instead of major updates being held back for the next release of the distribution. As these updates can bring major changes to one's system, it is possible for regular system updates to break a rolling release system.

The final category is semi-rolling. It is when some major components are allowed to roll, with major updates whenever they are ready, and other components stay at largely the same version (aside from security and bug fixes) until the next fixed release of the distribution. This approach comes with the pros and the cons of both fixed and rolling release distributions, with upgrades between fixed releases able to cause system breakage and so can regular updates at times. Although, regular updates, when they break things, usually are limited in the components they can break as not all major components are allowed to roll. 

# Alpine Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px" src="/Linux/Alpine_Linux_3.21.2.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">August 2005</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://alpinelinux.org/">alpinelinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">setup-alpine&mdash;text-based.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://wiki.alpinelinux.org/wiki/Alpine_Package_Keeper">Alpine Package Keeper</a> (APK; binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://wiki.alpinelinux.org/wiki/APKBUILD_Reference">APKBUILD</a>&mdash;shell script.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GNU Compiler Collection</a> (GCC)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/OpenRC">OpenRC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Musl">musl</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/BusyBox">BusyBox</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Almquist_shell">sh</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Detailed</td>
    </tr>
</table>
~~~

**Alpine Linux** is a security-focused distribution primarily intended for servers, routers, virtual private networks (VPNs), and alike. A base Alpine Linux install can be as small as 144 MB in size and does not include Bash, sudo, Vim or nano. The aforementioned intended uses are likely its ideal use case too, although I could see it being popular with desktop users that value security, a fast package manager, a minimalist system and a fixed release model. 

# Arch Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Arch_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">11 March 2002</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="http://www.archlinux.org/">www.archlinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer<exj/td>
        <td style="padding: 10px;"><a href="https://wiki.archlinux.org/title/Archinstall">archinstall</a>&mdash;textual installer.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://wiki.archlinux.org/title/Pacman">pacman</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://wiki.archlinux.org/title/PKGBUILD">PKGBUILD</a>&mdash;shell script.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a> (install), <a href="https://en.wikipedia.org/wiki/Z_shell">Zsh</a> (live).</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Vast, if the <a href="https://wiki.archlinux.org/title/Arch_User_Repository">Arch User Repository</a> (AUR) is included.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Comprehensive</td>
    </tr>
</table>
~~~
No conversation about Linux distributions geared towards advanced users would be complete without **Arch Linux**. It follows a Keep it Simple, Stupid (KISS) design philosophy. I may be biased in its favour as it is my go-to Linux distribution. A base install comes without a graphical user interface and has a pretty minimal array of packages, although the total size of a base install is about 1.7GB. It also has perhaps the most comprehensive documentation and vast repositories of any distribution, although it could be argued that NixOS has taken this title in recent years. That being said, I have experienced issues with Arch Linux before. Actually, I experienced one such issue while I was writing this webpage. See, I use Franklin.jl to build this website and I tried to deploy this website locally using my Arch Windows Subsystem for Linux (WSL) and I received errors related to the fact that Julia was using artefacts that expected OpenSSL 3.2.0 and my Arch WSL was using OpenSSL 3.4.0.

It is ideal for users that:

* Are comfortable with the command line. Those not comfortable with the command line may favour EndeavourOS or Manjaro Linux.
* Want the freedom to customize their system, but without the desire to compile most components of their system from source. 
* Want the very latest software. On the flip side of this, they should also know how to recover from an update breaking their system.
* Prefer a rolling release model.
* Prefer a fast package manager. pacman is one of the fastest I have ever encountered. 
* May want obscure pieces of software. Packaging on Arch is easy for people familiar with shell script &mdash; the language of the Linux command line &mdash; and with its vast repositories many users do not even need to resort to packaging software for themselves. 
* Do not mind using standard system software like systemd. Users that dislike systemd may prefer Artix Linux. 

# Chimera Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Chimera_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">2021</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://chimera-linux.org/">chimera-linux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Manual&mdash;bootstrapping.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://wiki.alpinelinux.org/wiki/Alpine_Package_Keeper">APK</a> (binary), <a href='https://github.com/chimera-linux/cports'>cports</a> (source)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://github.com/chimera-linux/cports">template.py</a>&mdash;Python script.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Clang">Clang</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://davmac.org/projects/dinit/">Dinit</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Musl">musl</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/FreeBSD">FreeBSD</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Almquist_shell">sh</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
</table>
~~~

**Chimera Linux** (not to be confused with ChimeraOS) is a truly unique Linux distribution and uses a very unusual combination of system software components. One interesting characteristic of the distribution that I did not mention in the infobox to the right is that Chimera Linux does not come with `sudo` pre-installed and it does not seem to be in Chimera's repositories ([source](https://pkgs.chimera-linux.org/packages?name=sudo&origin=)). Given the distribution's security focus, as evidenced by its use of musl, I would imagine this omission is a deliberate security feature. 

The ideal use case of Chimera Linux would be on security-critical systems, with users that favour FreeBSD's command line, do not need vast repositories and prefer rolling release models. Especially those that prefer to write their own packages using Python scripts, prefer fast package managers and dislike systemd.

# CRUX
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/CRUX_3.7.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">December 2002</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://crux.nu/">crux.nu</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Manual, with setup script.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">Ports with prt-get (source).</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://crux.nu/Main/Handbook3-7-Package">Pkgfile</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
</table>
~~~

**CRUX** aims to keep it simple as it uses tar.gz-based packages, BSD-style init scripts, and has fairly small repositories. It otherwise uses standard Linux system software. CRUX follows a fixed release model with new releases every year or two. It uses source-based package management and is best suited to advanced users that appreciate its idea of simplicity and want to compile their software from source. A base install of CRUX 3.7, with GRUB installed to serve as the bootloader, uses about 2.6GB disk space. 

# Debian
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Debian_12.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">August 1993</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://www.debian.org/">www.debian.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable<sup><a href="#footnote-5">5</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Debian-Installer">Debian-Installer</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">Advanced Packaging Tool</a> (APT; binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://wiki.debian.org/Packaging/Intro">Rules (Makefile), control, copyright and changelog files</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal or compete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Detailed</td>
    </tr>
</table>
~~~

**Debian** is the second-eldest Linux distribution still under active development. It has new stable releases every two years, roughly. It has three and sometimes four branches. In ascending order of modernity, they are (fourth branch in brackets): (old stable), stable, testing and unstable. Old stable corresponds to the previous stable release of the distribution. The stable branch corresponds to the current stable release of the distribution; each stable release comes with three years of support. In the lead up to a new stable release, the testing branch is forked and frozen and the packages undergo further testing and potentially patching until they are ready to be incorporated in the next stable release. Unstable is where Debian's very latest packages start out, until after sufficient testing they make their way into testing. Testing and unstable branches follow a rolling release model and are cutting edge and bleeding edge, respectively. 

Debian packages are built using a directory of packaging files. Among these is a rules file which is a Makefile with custom build commands. Personally, I have found Debian packaging one of the most challenging to wrap my head around. Partly because I found the custom build commands in rules files poorly documented. Although, naturally it may not be as challenging for others.

As users can choose a minimal install from its installer, and there are three main branches users can choose from, Debian can be a good choice for users that want to customize their system as much as one can without installing packages from source. Especially those that do not mind using systemd, like having very large repositories and do not mind having to use Makefiles to build packages, should one need to. Users needing more a beginner-friendly distribution should ideally use the Debian derivatives elementary OS, Linux Mint, MX Linux, Ubuntu or Zorin OS. 

Popular Debian derivatives include:
* [antiX](https://antixlinux.com/).
* [deepin](https://www.deepin.org/index/en) &mdash; although, in 2022 [they announced they were becoming independent](https://distrowatch.com/dwres.php?resource=showheadline&story=14870), but it seems to still be using Debian packages.
* [Devuan GNU+Linux](https://www.devuan.org/)
* [MX Linux](https://mxlinux.org/)
* [Ubuntu](https://ubuntu.com/) and its derivatives. 
I cover deepin and Ubuntu and some Ubuntu derivatives in separate sections, the rest I will cover here. 

antiX is designed to be lightweight and fast distribution with runit or SysV init as its init system. It uses JWM as its default user interface. It is ideal for users that want or need a lightweight distribution such as due to using old hardware. 

Devuan is essentially just Debian without systemd. It offers SysV init, runit and OpenRC editions.

# deepin
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/deepin_25Preview.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">28 February 2004</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://www.deepin.org/index/en">www.deepin.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Debian-Installer">deepin-Installer</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">APT</a> (binary) and <a href="https://www.deepin.org/en/deepin-linglong/">LingLong</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per Debian</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Compete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Detailed</td>
    </tr>
</table>
~~~
**deepin** uses its own desktop environment called [Deepin Desktop Environment](https://en.wikipedia.org/wiki/Deepin#Deepin_Desktop_Environment) (DDE) and has its own custom applications, including its own artificial intelligence (AI) assistant. It is developed by a Chinese company and has editions in both Mandarin Chinese and English (it seems to have support for other languages too, though). The English edition does have some untranslated Mandarin Chinese text in it, however. It is not enough to make the system unusable for people that do not understand this text, but may cause problems at times. Many consider it one of the most beautiful Linux distributions out there, at least in terms of its default aesthetics. Its packages can get outdated and it has tried to develop its own package manager called [Linglong](https://www.deepin.org/en/deepin-linglong/) as a way of providing more up to date versions of application software. 

It is ideal for users that want a beautiful desktop, have a large amount of free disk space, favour fixed release and appreciate distributions that try to innovate for their users. Especially users that are native Mandarin Chinese speakers. It also seems fairly beginner friendly to me. 

The free disk space required is at least 64 GB for the 25 preview release. The installer initially said 45 GB disk space was required, but when I went to partition my disk the installer said 64 GB disk space was required. The base installation ended up using just 18.6GB, roughly (after I had created my user account). I also noticed that deepin 25 preview used an immutable root file system. 

Its AI assistant answered my system memory when I asked how I could update my system given my root file system is read only. When I said that is not what I asked for it replied to me in Mandarin (even though my prompts were in English and my system language was set to English), although its answer this time seems relevant. Consequently, I would not say it is quite ready for everyday use, unless you have sometime to translate and correct its replies. 

# elementary OS
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/elementary_OS_8.0.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">31 March 2011</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://elementary.io/">elementary.io</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://github.com/elementary/installer">elementary Installer</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">APT</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per Debian.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**elementary OS** has its own desktop environment called [Pantheon](https://en.wikipedia.org/wiki/Elementary_OS#Pantheon_desktop_environment) which is built on [GTK 3](https://en.wikipedia.org/wiki/GTK) and [Vala](https://en.wikipedia.org/wiki/Vala_(programming_language)) and is rather aesthetically pleasing and has a macOS-like look with a dock. elementary OS is designed to be beginner-friendly and is based on Ubuntu LTS releases. macOS users that want to start try out Linux may prefer using elementary OS. It has a software centre that provides users the option to donate to their favourite projects. elementary OS itself can be downloaded for free, but its website does encourage users to pay what they want for the distribution.

# Exherbo
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Exherbo.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">17 January 2006?<sup><a href="#footnote-6">6</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://exherbolinux.org/">exherbolinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Manual&mdash;bootstrapping and compiling.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://paludis.exherbolinux.org/">Paludis</a> (source)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">ebuild&mdash;shell script with custom commands.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a><sup><a href="#footnote-7">7</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Exherbo** is a source-based distribution that originally was forked from Gentoo Linux. Like Gentoo, it uses ebuilds as its packaging files. Its package manager, Paludis, is written in [C++](https://en.wikipedia.org/wiki/C%2B%2B) unlike Gentoo's Portage, which is written in Python. Paludis is specifically meant to be a better alternative to Portage. Given Exherbo has smaller repositories and less comprehensive documentation, but is practically the same as Gentoo except without Gentoo's binary repositories, I would be inclined to think that Exherbo is best suited to Gentoo fans that are disgruntled with Portage. 

# Fedora
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Fedora_41.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">4 November 2003</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://fedoraproject.org/">fedoraproject.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Cutting edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Anaconda_(installer)">Anaconda</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/DNF_(software)">Dandified YUM</a> (DNF; binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Fedora** is a distribution that can be argued to be beginner friendly, although I am inclined to not put it in that category because it does not have out-of-the-box support for proprietary drivers, including WiFi drivers. Fedora is one of the most up-to-date fixed release distributions I am aware of, although each release usually keeps to the same release (except for bug and security fix releases) of desktop environment software and with six months between releases, this makes it not truly bleeding edge. Fedora releasers come with 13 months of support, so users only need to upgrade to every second release, should they choose. Fedora also has an immutable root file system edition called Silverblue. Fedora is best suited to users that favour a fixed release model, like cutting edge software, need large repositories and prefer to package with spec files, when this is necessary. 

Fedora is the basis of [CentOS Stream](https://en.wikipedia.org/wiki/CentOS_Stream), which in turn is the basis of [Red Hat Enterprise Linux](https://en.wikipedia.org/wiki/Red_Hat_Enterprise_Linux) (RHEL) and derivatives thereof like [AlmaLinux](https://en.wikipedia.org/wiki/AlmaLinux) and [Rocky Linux](https://en.wikipedia.org/wiki/Rocky_Linux). RHEL and its derivatives are popular server distributions; they come with about a decade of support. RHEL itself comes with an additional two years of extended lifecycle support. 

# Gentoo Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Gentoo_Linux.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">31 March 2002</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://www.gentoo.org/">www.gentoo.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Cutting edge<sup><a href="#footnote-8">8</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Manual&mdash;bootstrapping and compiling.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Portage_(software)">Portage</a> (source)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">ebuild&mdash;shell script.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/OpenRC">OpenRC</a>/<a href="https://en.wikipedia.org/wiki/Systemd">systemd</a><sup><a href="#footnote-7">7</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Detailed</td>
    </tr>
</table>
~~~

**Gentoo Linux** is, in many ways, the prototypical source-based Linux distribution. It uses its own package manager called Portage which is meant to be a true ports system in the spirit of BSD ports. Although, in recent years it has become feasible to install most package as pre-compiled binaries via Portage. Interestingly, there have even been projects to port Portage to other operating systems like the BSD derivatives FreeBSD and NetBSD.

Previously, [Sabayon Linux](https://en.wikipedia.org/wiki/Sabayon_Linux) occupied this niche by offering binary packages on a Gentoo base while still allowing users to install software from source via Portage. Although, Sabayon provided binary packages via its own package manager called Entropy. Sabayon Linux was discontinued around 2020. 

[Calculate Linux](https://en.wikipedia.org/wiki/Calculate_Linux), which is still actively maintained, can also be argued to occupy this niche as it provides binary packages too while still giving users the option to install from source via Portage. Unlike Sabayon, which provided its own binary package manager, Calculate just uses Portage to install binary packages. Sabayon and Calculate both have or had automated installers, unlike Gentoo. 

Now it seems like Gentoo itself is trying to occupy the niche of offering binary packages on a Gentoo base as well. As a casual user that likes to try it out in virtual machines from time to time, is an attractive feature in my opinion.

Gentoo is ideal for Linux users that want complete freedom to customize their system all the way down to the configure/compile options used to build each package. Users can even fork packages and apply custom patches to them, should they choose. 

# Guix System
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Guix_System.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">29 March 2016<sup><a href="#footnote-9">9</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://guix.gnu.org/">guix.gnu.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Text-based installer.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Guix">GNU Guix</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Guile">GNU Guile</a> scripts.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Guix#GNU_Shepherd_init_system">GNU Shepherd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Detailed</td>
    </tr>
</table>
~~~

**Guix System** (pronounced "Geeks") is a reproducible, entirely free (as in freedom) Linux distribution that uses its own package manager called GNU Guix which installs each package to its own unique directory within `/gnu/store`. While Guix System uses GNU Guix as its package manager, GNU Guix is technically distribution-agnostic. Guix System is configured using files written in GNU Guile, such as `/etc/config.scm`. GNU Guile is also used to write packaging files for GNU Guix. Unlike NixOS, another reproducible Linux distribution, it does not seem to keep old configurations in its bootloader menu. It seems suitable for users that want a system entirely configurable using a single file written in GNU Guile and favour an entirely free operating system, even though this often comes with hardware compatibility issues. 

# Linux From Scratch
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Linux_From_Scratch.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">3 December 1999</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="http://www.linuxfromscratch.org/">www.linuxfromscratch.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Manual compilation of each component.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">None, software manually compiled from source.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">None</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a>/<a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Detailed</td>
    </tr>
</table>
~~~

**Linux From Scratch** (LFS) is a source-based distribution wherein each software package is manually compiled and installed from source. Users achieve this by following the instructions in a book provided by the LFS project. LFS itself only provides users with a base Linux system, there is a sister project called [Beyond Linux From Scratch](https://www.linuxfromscratch.org/blfs/) (BLFS) that provides users with the additional software (e.g. graphical user interface software) required for a more complete and functional system. LFS does not have a piece of software to manage package management for the user, instead the user is the package manager. This does give users complete ability to build their system from the ground up and customize it to their liking. 

Many people find installing LFS a frustrating experience as it is tedious and small errors can cause big problems. Despite these frustratons, installing a LFS system is a very effective way to learn about the inner workings of a Linux operating system. It is also an invaluable option, as far as Linux distributions go, for Linux users that want to customize their system down to the compile options from the ground up. 

# Linux Mint
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Linux_Mint_22.1.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">27 August 2006</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://linuxmint.com/">linuxmint.com</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Ubiquity&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">APT</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per <a href="#debian">Debian</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Linux Mint** is a beginner-friendly distribution based on Ubuntu's long-term support (LTS) releases. Its team forked GNOME 3 to create Cinnamon in an attempt to provide users a more classic desktop experience. It has three official editions that all feature a classic desktop experience that includes a Windows-like layout. The distribution includes many of its own tools for common tasks like package management and configuration. There is also a Debian-based edition of Mint. 

It is ideal for beginners that are used to the Windows layout, especially if they do not want the latest software, would rather have system upgrades every two years or so, and do not have especially exotic software needs. 

# Mageia
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Mageia_9.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">1 July 2011</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://www.mageia.org/en">www.mageia.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">DrakX&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">DNF (current) and urpmi (legacy)&mdash;both binary.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Spec file.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Mageia** is a Linux distribution that started out in 2011 as a fork of Mandriva Linux created by some former employees of the company that had developed Mandriva. Originally, it used the same tools as Mandriva like the package manager urpmi, but it has modernized in some ways and now uses DNF as its package manager. My experience with it is that it is rock solid stable, but many packages that I use are missing from their repositories. Consequently, I would recommend Mageia to users that want a rock solid stable system and do not have obscure software needs. Especially if they used Mandriva Linux and were fond of it.  

# MX Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/MX_Linux_23.5.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">24 March 2014</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://mxlinux.org/">mxlinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">MX Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">Advanced Packaging Tool</a> (APT; binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per Debian</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Compete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~
**MX Linux** is based on antiX but uses customized and prettified Xfce, KDE Plasma or Fluxbox as its user interface. It has several tools specifically developed for the distribution, including configuration tools and graphical package management tools. It is ideal for users that dislike systemd, like aesthetically pleasing desktops, and prefer a fairly beginner-friendly experience. 

# NixOS
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/NixOS_24.11.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">3 June 2003</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://nixos.org/">nixos.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Nix_(package_manager)">Nix</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Nix file</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Vast</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Comprehensive</td>
    </tr>
</table>
~~~

**NixOS** is a reproducible Linux distribution that is configured using files written in a special-made purely functional programming language called [Nix](https://nix.dev/manual/nix/latest/language/index.html). While it is purely functional, it does incorporate a few syntactic elements from shell script. It also uses a package manager called [Nix](https://nix.dev/manual/nix/2.25/introduction); Nix installs each package to its own unique directory within `/nix/store`. This means that multiple versions of the same package can be installed on a NixOS system. It then sets up symlinks and environment variables to ensure that each piece of software is able to find all libraries, binaries and alike that it depends on. Nix packages are also specified using files written in the Nix language. Nix and NixOS started out in the early 2000s as a research project by then software engineering student Eelco Dolsta. 

Its chief system configuration file is [`/etc/nixos/configuration.nix`](https://nixos.org/manual/nixos/stable/#ch-configuration) and this file largely uniquely determines the root file system of the distribution. This is why the system is reproducible, as the root file system of two NixOS systems built with the same configuration file will be largely the same. This is with the exception that if additional packages are installed using user configuration files or running `nix-env -i <package>` they will be installed under `/nix/store`. Whenever one wants changes to the aforementioned system configuration file to come into effect, one runs `nixos-rebuild switch` (as root) and the new configuration is built. The old configuration is also kept and when users boot the system they can boot the new configuration (which is the default), or the old configuration. NixOS also keeps even older configurations, if they exist, although naturally this uses disk space so there is a command to remove older configurations (`nix-collect-garbage -d`) to free up disk space. 

One thing I like about NixOS is that it usually will not let you build an invalid configuration, which means that whenever I boot NixOS, I can be almost certain it will successfully boot. NixOS has one other system configuration file that I have not mentioned, although it specifically pertains to hardware configuration. It is `/etc/nixos/hardware-configuration.nix` and it is where I have found NixOS seems to turn a blind eye to certain errors and sometimes will let you build an invalid configuration. Specifically, I have built unbootable NixOS systems by accident by specifying a root file system in this file that does not exist. No warning was given that I had specified an invalid root file system. I have started a discussion about this issue on [NixOS' discourse](https://discourse.nixos.org/t/automatic-hardware-detection-for-nixos-rebuild-nixos-install/59666).

NixOS is ideal for intermediate to advanced Linux users that:

* Want a reproducible system.
* Do not mind using systemd.
* Like the idea of configuring their system using a file written in a functional programming language.
* Want a system that it is more difficult to break, as it typically will not allow you to build a broken configuration. 
* Would like to create packages for their system using files written in that language.

Nix has also been ported to several other operating systems, including BSD derivatives and macOS. There is even a [NixBSD project](https://github.com/nixos-bsd/nixbsd) that aims to create an operating system that uses the Nix package manager and the FreeBSD kernel. 

# NuTyX
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/NuTyX_24.10.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">14 September 2009<sup><a href="#footnote-10">10</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://nutyx.org/en/">nutyx.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">NuTyX Installer&mdash;text-based.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">Cards (binary and source)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Pkgfile&mdash;shell script</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
</table>
~~~

**NuTyX** is based on LFS and BLFS but has its own package manager called cards. It allows users to install software from binary packages and from source via a ports system. This makes it remind me a little of FreeBSD's approach to package management, as it has a binary package manager called pkg and a ports system that users can use to install software from source. NuTyX is aimed at intermediate to advanced users. I think NuTyX is ideal for intermediate to advanced users that do not need obscure software, and want a distribution with a hybrid approach to package management.

Something interesting about NuTyX that I noticed in one of my virtual machines that runs NuTyX 24.10 was that installing VirtualBox guest additions seems to cause a few errors to appear in the boot screen, although the system still successfully boots and loads a graphical user interface. 

# openmamba GNU/Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/openmamba_GNU_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">Before 4 August 2009<sup><a href="#footnote-11">11</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://openmamba.org/">openmamba.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Cutting edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Mamba Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">DNF (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
</table>
~~~

**openmamba GNU/Linux** is a Linux distribution that offers out-of-the-box support for hardware with proprietary drivers. I personally found it ran fine in a virtual machine, but I have seen some reviews of it online that have mentioned significant bugs in previous installation medium releases. I am inclined to suggest it as an option for users fond of RPM packaging and rolling release models that do would be content with the distribution's relatively small repositories and need out-of-the-box support for devices that require proprietary drivers.

# OpenMandriva Lx
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/OpenMandriva_Lx_Rolling.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">18 June 2013</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://www.openmandriva.org/">www.openmandriva.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed and rolling</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable (fixed), bleeding edge (rolling)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">urpmi (legacy) and DNF (current)&mdash;both binary.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Clang">Clang</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**OpenMandriva Lx** is a continuation of Mandriva Linux developed by a community project. Like Mageia, it uses the DNF package manager. Unlike Mageia, it comes in two editions&mdash;a fixed release and rolling release edition. It first started to offer a rolling release edition in 2023.~~~<sup><a href="#footnote-12">12</a></sup>~~~ One major difference with Mageia is that it uses Clang as its compiler. It seems most suitable to users that favour RPM packaging, want to use a distribution with Clang-compiled packages, do not mind its relatively small repositories and have a fondness for the old Mandriva Linux distribution.

# openSUSE
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/openSUSE_Tumbleweed.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">7 December 2006</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://www.opensuse.org/">www.opensuse.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed (Leap) and rolling (Tumbleweed)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable (Leap), bleeding edge (Tumbleweed)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">YaST&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">ZYpp (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**openSUSE** is a continuation of the SUSE Linux distribution developed by a group of German computer science students and first released in March 1994. Like OpenMandriva Lx, it comes into separate editions &mdash; one, Leap, that features a fixed release model and another, Tumbleweed, that features a rolling release model. openSUSE started providing two separate editions in 2014,~~~<sup><a href="#footnote-13">13</a></sup>~~~ whereas OpenMandriva Lx adopted this two edition approach around 2023.~~~<sup><a href="#footnote-12">12</a></sup>~~~ 

One notable feature of openSUSE is that, by default, it uses [Btrfs](https://en.wikipedia.org/wiki/Btrfs) as its root file system. It is used as it allows for easier snapshots to backup the root file system. In my experience, this is more of a curse than a blessing, as I tend to find that openSUSE with a Btrfs root file system becomes unbootable within about a fortnight for me, at least. This is even when I keep on top of the snapshots, delete the old ones and keep an eye on my disk usage using Btrfs' own tools. I mention that I use Btrfs' own tools as the Linux command-line tool `df` is not accurate in measuring file system usage when it is a Btrfs file system. openSUSE also uses [XFS](https://en.wikipedia.org/wiki/XFS) as its default home file system. 

I would recommend openSUSE to intermediate to advanced users that like RPM packaging and may need obscure pieces of software. In theory, it could be used by a beginner, but I personally think that a beginner would likely really struggle with Btrfs. This is obviously a problem for openSUSE given it is the default root file system of the distribution. Either they would struggle to keep on top of the snapshots and preventing them from using up their entire root file system, or they may experience an issue like that one I previously mentioned. 

[SUSE Linux Enterprise](https://en.wikipedia.org/wiki/SUSE_Linux_Enterprise) (SLE) is based on openSUSE Leap and is a commercial product. It comes with about thirteen years of general support, an additional three years of long term service pack support (total sixteen years) and an additional three years of long term service pack core support (total nineteen years). It comes in two editions, SLE Server for servers and mainframes and a desktop/workframe edition called SLE Desktop.

# PCLinuxOS
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/PCLinuxOS.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">October 2003</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://pclinuxos.com/">pclinuxos.com</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">PCLinuxOS Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/APT-RPM">APT-RPM</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Spec file</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**PCLinuxOS** is a beginner-friendly Linux distribution that was originally forked by Bill Reynolds (Texstar) from Mandrake Linux 9.2 in 2003. It is rather conservative in some ways, for instance it still uses SysV as its init system, APT-RPM as its command-line package manager and Synaptic as its graphical package manager. APT-RPM had its last release in 2008, Synaptic has an outdated look although it is still maintained and SysV has been superseded on most distributions (not all, of course) by systemd, which was first released in 2010. Despite using a rolling release model, it also usually uses pretty old software. 

PCLinuxOS is perhaps best suited to beginners that do not need exotic software, like a no thrills and 2000s-esque desktop experience and favour a rolling release model. If somehow, despite being beginners, they have an opinion on init systems and dislike systemd, they may also like PCLinuxOS. 

# Pop!_OS
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><image width="380px;" src="/Linux/Pop!_OS_22.04.png"/></td>
    </tr> 
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">27 October 2017</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://pop.system76.com/">pop.system76.com</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Old stable<sup><a href="#footnote-14">14</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://github.com/elementary/installer">elementary Installer</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">APT</a> (binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per Debian</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Compete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">None<sup><a href="#footnote-15">15</a></sup></td>
    </tr>
</table>
~~~

**Pop!_OS** originally used a customized GNOME desktop but its team has been developing a desktop environment written in Rust called COSMIC. It is beautiful by default. Pop!\_OS is developed by the computer manufacturer called System76. According to its website, it is aimed at STEM and creative professionals. It does seem fairly beginner friendly from my experience with it, but it does use pretty old software due to it, at the time of writing (26 January 2025), being based on the previous long-term support (LTS) release of Ubuntu. I would recommend Pop!\_OS to users that want an eye candy desktop by default, do not mind older packages and prefer a fixed release distribution.

# Rhino Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Rhino_Linux.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">8 August 2023</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://rhinolinux.org/">rhinolinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://github.com/volitank/nala">Nala</a> (binary), <a href="https://github.com/pacstall/pacstall">Pacstall</a> (source) and distro-agnostic package managers (binary).<sup><a href="#footnote-16">16</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per Debian, plus Pacscript (shell script).</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Vast</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Rhino Linux** is a Ubuntu derivative that gets its core packages from the Ubuntu developmental branch. Rhino Linux is the only Ubuntu-based distribution that I am aware of that follows a rolling release model. It is specifically designed with developers in mind and comes with VSCodium pre-installed. 

It also comes with pacstall &mdash; a package manager that provides access to a repository designed to be the Ubuntu counterpart to the Arch User Repository &mdash; pre-installed. Rhino also has a setup wizard that offers users four different distribution-agnostic package formats that the wizard can add support for onto their system &mdash; namely Snap, Flatpak, Nix and AppImages (with the AppImage manager [AM](https://github.com/ivan-hc/AM)). Rhino has a wrapper for each of its package managers (including the distribution-agnostic package managers) that is called rhino-pkg. 

Its default desktop environment is a customized Xfce desktop featuring a dock on the left of the screen along with a global menu. This customized desktop they call "Unicorn". Unicorn has a default look that features a lot of purple, black and white that appeals to my eye, at least. 

I would say that Rhino Linux is probably ideal for developers. Especially those that:

* Prefer a rolling release model.
* Prefer graphical approaches to installation and package management.
* Have obscure software needs.
* Like eye candy distributions.
* Prefer the Xfce desktop.
* Like Ubuntu-based distributions.
* Prefer using shell script to package for their system.
* Want the very latest software. 

I personally rather like this distribution as it addresses many of the problems I had with Ubuntu when I stopped using it as my daily driver around 2015. 

# Slackware Linux
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Slackware_Linux_15.0.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">17 July 1993</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="http://www.slackware.com/">www.slackware.com</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Slackware Installer&mdash;Textual.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;">pkgtools</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://docs.slackware.com/slackware:slackbuild_scripts">SlackBuilds</a>&mdash;shell script</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Slackware Linux** is the eldest Linux distribution still in active development. An important characteristic of it is that its official repositories are fairly small and mostly just contain the packages that one can install from the live medium and updates thereof. There are unofficial repositories, but even they are not very large. This is largely because on Slackware it is expected that most non-core packages will be manually compiled from source using SlackBuild scripts. 

Another important characteristic of Slackware is that its developers are fairly conservative in that they are reluctant to adopt divisive pieces of technology like systemd, and often ship pretty old and well-tested versions of the software included in the system. I say they, but Slackware technically has a Benevolent Dictator for Life named Patrick Volkerding, who was its original creator back in 1993. It is also worth noting that stable releases of Slackware have become fairly rare. The period between release 14.2 and 15.0, for instance, was 5 years and 7 months, roughly. 

I have tried Slackware many times and I have found its approach to package management frustrating. That being said, it is rock solid stable and if you are nostalgic for how Linux distros were like in the 1990s and like its approach to package management, it may be a suitable distribution for you. Especially if you are dislike the inclusion of systemd in most other modern Linux distributions.

# Solus
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Solus.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">27 December 2016</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://getsol.us/">getsol.us</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Calamares_(software)">Calamares</a>&mdash;Textual.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Pardus_(operating_system)#PiSi_package_management">eopkg</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://help.getsol.us/docs/packaging/package.yml/">package.yml</a>&mdash;YAML file.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Init#SysV-style">SysV</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal or complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

Like PCLinuxOS, **Solus** is a Linux distribution aimed towards beginners, despite following a rolling release model. As its package manager, it uses eopkg&mdash;which is based on [Pardus](https://en.wikipedia.org/wiki/Pardus_(operating_system))' abandoned package manager of PiSi. It is the original distribution that the Budgie desktop was developed for, and is noted for its relatively good default aesthetics. 

Its ideal use case is probably a beginner that does not need exotic software, appreciates a beautiful and simple desktop like Budgie and does not want to have to upgrade their system between releases of a fixed release distribution like Linux Mint for fear of system breakage. Granted, system updates on a rolling release system can break things too, so this should be factored in when considering PCLinuxOS. 

# Ubuntu
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Ubuntu_24.10.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">20 October 2004</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://ubuntu.com/">ubuntu.com</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Ubiquity_(software)">Ubiquity</a>&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Advanced_Packaging_Tool">APT</a> (binary) and <a href="https://en.wikipedia.org/wiki/Snap_(software)">Snap</a> (binary).</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Per <a href="#debian">Debian</a>.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a><td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Very large</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

**Ubuntu** is a Linux distribution based on Debian (testing) and has been the go-to beginner-friendly distribution since its first release in 2004. It was created by South African enterpeneur Mark Shuttleworth and is maintained by his company [Canonical](https://en.wikipedia.org/wiki/Canonical_(company)). Many consider Ubuntu responsible for a lot of the changes in the Linux world that have made Linux distributions more accessible to novice users. Canonical has also been an innovator in other ways, such as by developing the Ubiquity system installer, [Upstart init system](https://en.wikipedia.org/wiki/Upstart_(software)), the Snap distribution-agnostic package manager, [Mir display server](https://en.wikipedia.org/wiki/Mir_(software)) and [Unity graphical shell](https://en.wikipedia.org/wiki/Unity_(user_interface)). Although, Upstart is no longer developed and Unity is no longer under development by Canonical. Ubuntu was an early adopter of each of these technologies.

New Ubuntu releases come out every six months, usually in April and October of every year since its initial release in October 2004. In April of even-numbered years, there are long-term support (LTS) of the distribution that receive about five years of support. Other releases receive nine months of support. The nine monthly releases usually come with the latest desktop environment releases or near to it, and a fairly modern kernel. 

Ubuntu is ideal for beginners that favour a fixed release cycle. Given its two types of editions one with long support periods, it gives users a lot of choice for when they will need to upgrade their system.

Ubuntu is a very popular base for other distributions. I am not going to cover all Ubuntu derivatives, some that I will not cover separately that are of note are [KDE neon](https://neon.kde.org/), [Linux Lite](https://www.linuxliteos.com/), [TUXEDO OS](https://www.tuxedocomputers.com/os) and [Zorin OS](https://zorin.com/os/). I am not covering these distributions separately as I do not see them as innovative enough to warrant it. This is not to insult the developers, they are perfectly acceptable distributions to use and I can definitely see work that went into them, but I do not have enough to mention about them to warrant a separate section. 

KDE neon is a semi-rolling release distribution as its core system software is based on Ubuntu LTS releases but its KDE software is bleeding edge. It is not especially beginner friendly and its ideal users are KDE fans that want to try out the latest KDE software as soon as it is published.

Linux Lite is a beginner-friendly and fairly light Ubuntu LTS-based distribution. It has a beautiful default look that resembles Windows. 

TUXEDO OS is a beginner-friendly Ubuntu LTS-based distribution developed by TUXEDO Computers in Germany. It has its own control centre.

Zorin OS uses a Windows-like layout and is aimed at beginners as well. 

# Vanilla OS
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Vanilla_OS_2.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">29 December 2022<sup><a href="#footnote-17">17</a></sup></td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://vanillaos.org/">vanillaos.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;">Fixed</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Stable</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Vanilla Installer&mdash;graphical.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/OSTree">OSTree</a> (read-only root), distro-agnostic package managers<sup><a href="#footnote-18">18</a></sup> and <a href="https://github.com/Vanilla-OS/apx">Apx</a> (applications).</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;">Any, due to Apx.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Systemd">systemd</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a><td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Vast</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Complete</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Medium</td>
    </tr>
</table>
~~~

I would like to start this section by mentioning that I have not been able to install **Vanilla OS**, as the installer repeatedly fails for me. That being said, Vanilla OS uses a Debian (unstable) base (although, it previously used a Ubuntu base) and an immutable root file system. Actually, it has two root file systems. One is booted by the user, the other is the one to which updates are applied. This is so that users can boot this updated system at their next reboot but have their other root file system as a backup should the update break their system. Vanilla OS uses Apx (pronounced "apex") which uses containerized Linux distributions to provide access to software packaged for that distribution. 

Vanilla OS sounds ideal for at least intermediately experienced users that have plenty of available disk space, want a Debian-based immutable system and access to packages that Debian does not provide. The requirement for a lot of disk space is due to the two root file systems and containerized Linux distributions to provide additional packages, both of which substantially more space than most distributions would need.

# Void
~~~
<table style="width: 380px; float: right; border-collapse: collapse;">
    <tr>
        <td colspan="2"><img width="380px;" src="/Linux/Void.png"/></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Initial release</td>
        <td style="padding: 10px;">2008</td>
    </tr>
    <tr>
        <td style="padding: 10px; width: 190px;">Website</td>
        <td style="padding: 10px;"><a href="https://voidlinux.org/">voidlinux.org</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Release model</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Rolling_release">Rolling</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Modernity<sup><a href="#footnote-1">1</a></sup></td>
        <td style="padding: 10px;">Bleeding edge</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Installer</td>
        <td style="padding: 10px;">Void Installer&mdash;textual.</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Package manager (type)</td>
        <td style="padding: 10px;"><a href="https://github.com/void-linux/xbps">X Binary Package System</a> (XBPS; binary)</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Packaging file(s)</td>
        <td style="padding: 10px;"><a href="https://xbps-src-tutorials.github.io/">template</a>&mdash;shell script</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Compiler</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU_Compiler_Collection">GCC</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Init system</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Runit">runit</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">C standard library</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Glibc">glibc</a>/<a href="https://en.wikipedia.org/wiki/Musl">musl</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Userland</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/GNU">GNU</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Shell</td>
        <td style="padding: 10px;"><a href="https://en.wikipedia.org/wiki/Bash_(Unix_shell)">Bash</a></td>
    </tr>
    <tr>
        <td style="padding: 10px;">Repository size<sup><a href="#footnote-2">2</a></sup></td>
        <td style="padding: 10px;">Medium-small</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Base install<sup><a href="#footnote-3">3</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
    <tr>
        <td style="padding: 10px;">Documentation<sup><a href="#footnote-4">4</a></sup></td>
        <td style="padding: 10px;">Minimal</td>
    </tr>
</table>
~~~

**Void** is a Linux distribution that is similar to Arch Linux in that it follows a rolling release model and uses a fast, lightweight package manager written in C that uses shell script packaging files for building its packages. A base install without a GUI is about 3GB in size. It has a fairly small development team, but I find it an interesting system. It boots rather fast by default and it has a nice air of Arch Linux-style simplicity. 

Void is ideal for experienced users that are content with the software in Void's relatively small repositories, prefer shell script for packaging, want a fast package manager, do not mind a command-line installation process, prefer runit to systemd, and would prefer the option to use musl instead of glibc. 

# Footnotes
~~~
<ol>
    <li id="footnote-1">The possible categories are, in ascending order of modernity: old stable, stable, cutting edge and bleeding edge.</li>
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
    <li id="footnote-14">This is based on the fact that the latest release as of 26 January 2025 is based on Ubuntu 22.04.</li>
    <li id="footnote-15">I checked the Pop!_OS website and could not find documentation on it.</li>
    <li id="footnote-16">I mention these, even though I omit them in most distro's infoboxes, because Rhino Linux has options to enable cross-distro package managers/formats in its setup wizard. Specifically, it allows users to enable Flatpak, Nix, Snap or AppImages.</li> 
    <li id="footnote-17"><a href="https://vanillaos.org/blog/article/2022-12-29/vanilla-os-2210-kinetic-is-out">Vanilla OS 22.10</a> was this release and it was the first release mentioned in the <a href="https://vanillaos.org/blog">Vanilla OS</a> blog.</li>
    <li id="footnote-18">I say this because <a ref="https://docs.vanillaos.org/handbook/en/install-flatpaks">Flatpak installation instructions</a> and <a href="https://docs.vanillaos.org/handbook/en/install-homebrew">Homebrew installation instructions</a> are given in Vanilla OS's documentation.</li>
</ol>
~~~