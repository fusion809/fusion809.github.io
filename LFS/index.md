@def hassim=false;
@def title="Linux From Scratch (LFS) virtual machine (VM)"
@def tag=["Linux"]
@def mintoclevel=1

![LFS screenshot](https://fusion809.github.io/images/executor-raujonas.github.io/LFS_screenshot_02-10-2026.png)

**Figure 1: Screenshot of my LFS VM's GNOME session as of 2 October 2026.**

I first installed LFS 12.4 systemd edition to a virtual machine on 9 February 2026. Since then, I have upgraded the system to the development systemd branch, and then gradually made it even more bleeding edge than this by upgrading all packages to the latest stable upstream release. Sometimes I need to keep a package back simply because its latest stable release actually depends on pre-release versions of other packages. This was the case for `gnome-control-center` on 12 September, as at this point version `51.0` of the package was out but it depended on version `51.alpha` or later of `gnome-desktop` and out of these versions only `51.alpha` was available at the time. 

\toc

\includemd{LFS/_motivations.md}

\includemd{LFS/_package-management.md}

\includemd{LFS/repos/_overview.md}

\includemd{LFS/repos/_build_duration.md}

\includemd{LFS/repos/_apps.md}

\includemd{LFS/repos/_dotfiles.md}

\includemd{LFS/repos/_logs.md}

\includemd{LFS/repos/_packaging.md}

\includemd{LFS/repos/_scripts.md}

\includemd{LFS/_gnome.md}

\includemd{LFS/_kde.md}
