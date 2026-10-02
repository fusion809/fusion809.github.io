# Package management
From my NixOS host machine, I have written &mdash; with the help of artificial intelligence (AI) &mdash; several shell functions that are imported into my LFS VM and provide basic package management functionality. These functions are part of both my host's and VM's shell profile. These functions can be found in my [NixOS configuration user shell profile](https://github.com/fusion809/NixOS-configs/tree/26.05/shell/user/). [lfs-custom-updates.py](https://github.com/fusion809/NixOS-configs/blob/26.05/python/lfs-custom-updates.py) is used to parallelize and more efficiently check the versions of all custom packages to see if updates are available.

~~~
<table style="border-collapse: collapse; width: 100%;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 1: Shell functions used for package management within the LFS VM.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center; white-space: nowrap;">
            <b>Syntax</b> (definition file hyperlink)
        </td>
        <td style="font-size: 20px; padding: 10px; overflow-wrap: break-word;">
            <b>Description</b>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autobuild-func.sh"><code>autobuild PACKAGE(S) [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Default: build and install the specified package(s), if and only if the latest version of the package is not already installed. LFS/Beyond LFS (BLFS) instructions were used to build most packages. Although, some packages were built using custom build scripts defined in <a href="https://github.com/fusion809/lfs_packaging"><code>~/lfs_packaging</code></a>. These custom build scripts have since become the primary source of software for the system, as it tends to be more reliable and easier to tweak.<br/>
            Options:<br/>
            <code>--dry-run</code>: show what actions would be executed to build and install the package.<br/>
            <code>--strip</code>: run stripping commands after build.<br/>
            <code>--no-upstream</code>: disable upstream version searching.<br/>
            <code>--include-config</code>: include configuration commands in the LFS/BLFS book entry.<br/>
            <code>--rm-libs</code>: remove old library versions after build (disabled by default).<br/>
            <code>--lfs</code>: search only in the LFS book.<br/>
            <code>--blfs</code>: search only in the BLFS book.<br/>
            <code>--lfs-book &lt;book&gt;</code>: specify LFS book (e.g., development, systemd, stable, or full URL).<br/>
            <code>--blfs-book &lt;book&gt;</code>: specify BLFS book (e.g., development, systemd, stable, or full URL).<br/>
            <code>--skip-tests</code>: skip test commands (make check/test, etc.).<br/>
            <code>--ignore-test-failures</code>: ignore test failures by appending '|| true' to test commands.<br/>
            <code>-f</code>/<code>--force</code>: force rebuild and installation even when latest version is already installed.<br/>
            <code>-h</code>/<code>--help</code>: show help message.<br/>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-autoremove.sh"><code>autoremove PACKAGE(S) [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Default: remove the specified package(s), if and only if no other packages have libraries that depend on the package(s).<br/>
            Options:<br/>
            <code>--dry-run</code>: show what actions would be executed to remove the package.<br/>
            <code>-f</code>/<code>--force</code>: force removal, without regard for library dependencies.<br/>
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-libs.sh"><code>ls_old_libs [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            List old installed versions of libraries.<br/>
            Options:<br/><code>-d</code> option it list files that depend on listed files.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/Shell/02-pms.sh"><code>rm_book_src</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove book source files.
        </td>
    </tr>  
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/Shell/02-pms.sh"><code>rm_lfp_src</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove custom package source tarballs.
        </td>
    </tr>    
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-share.sh"><code>rm_old_docs</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old unused documentation directories.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-kerns.sh"><code>rm_old_kerns</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old unused kernels.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-libs.sh"><code>rm_old_libs</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old unused libraries. As for used old used libraries, rebuild packages that depend on the library and then remove it.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-share.sh"><code>rm_old_share [OPTION(S)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Removes old unused <code>/usr/share</code> subdirectories.<br/>
            Options:<br/><code>--dry-run</code> shows what would be done without actually executing those actions.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/Shell/02-pms.sh"><code>rm_src</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Remove old source archives and directories (not including git repos).
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-sync-to-vm.sh"><code>sync_to_vm</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Only available on host; synchronize scripts from host to virtual machine. 
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-update.sh"><code>update [OPTION(s)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Update packages.<br/>
            Options:<br/>
            <code>--dry-run</code>: Show what would be updated without downloading/building.<br/>
            <code>-h</code>/<code>--help</code>: Show help message.<br/>
            <code>--no-upstream</code>: Check only LFS/BLFS book versions (disable upstream tracking).<br/>
            <code>-v</code>/<code>--verbose</code>: Show custom package local and upstream versions during version checking.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-update.sh"><code>updatec [OPTION(s)]</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Runs, in order, <code>rm_old_docs</code>, <code>rm_old_kerns</code>, <code>rm_old_libs</code>, <code>rm_old_share</code> and <code>rm_src</code> if and only if <code>update</code> runs without error. Options are passed directory to <code>update</code>.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px; white-space: nowrap;">
            <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh"><code>updates</code></a>
        </td>
        <td style="font-size: 16px; padding: 10px; overflow-wrap: break-word;">
            Print table that list available updates (marked with <code>[UPDATE]</code>), as well as packages with missing inventories (marked with <code>[FILES MISSING]</code>) and packages with versioning failures (marked with <code>[FAILED]</code>). 
        </td>
    </tr>
</table>
~~~
