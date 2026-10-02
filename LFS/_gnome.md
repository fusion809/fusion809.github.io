# GNOME
GNOME was the first desktop I installed and is the main user interface I use in the virtual machine. My NixOS system uses Hyprland instead, but I have struggled to get Hyprland to actually work in a KVM/QEMU virtual machine, so I decided to just use GNOME in the LFS VM. 

[Dash to Dock](https://github.com/micheleg/dash-to-dock) is enabled and installed, as is [WeatherPanel](https://github.com/attentivecoder/weatherpanel), [Extension List](https://github.com/tuberry/extension-list), [Kiwimenu](https://github.com/kem-a/kiwi-menu), [Show Desktop Button](https://github.com/amivaleo/Show-Desktop-Button) and [Super Into Apps](https://github.com/mikelei8291/super-into-apps). As previously mentioned, I also use my own [own fork](https://github.com/fusion809/executor-raujonas.github.io) of the Executor extension. 

~~~
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 2: GNOME themes.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Cursor</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Legacy applications</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Icons</b>
        </td>
        <td style="font-size: 20px; padding: 10px;">
            <b>Shell</b>
        </td>
    </tr>
    <tr>
    <td style="font-size: 16px; padding: 10px;">
    WhiteSur-cursors
    </td>
    <td style="font-size: 16px; padding: 10px;">WhiteSur-Dark
    </td>
    <td style="font-size: 16px; padding: 10px;">WhiteSur-dark
    </td>
    <td style="font-size: 16px; padding: 10px;">TST - Semi Transparent
    </td>
    </tr>
</table>
~~~

## Executor fork
*Has `~/lfs_packaging` package called [executor](https://github.com/fusion809/lfs_packaging/tree/master/executor). It can be installed via more standard ways, too.*

The base [Executor](https://github.com/raujonas/executor) extension provides up to three widgets in the GNOME panel on the left, centre and right of the panel. In these widgets is displayed the output of specified commands. The interval at which the command is re-run can also be specified. The [Executor fork](https://github.com/fusion809/executor-raujonas.github.io) I maintain provides the following additional features:
* Tooltips &mdash; which can have two components. They are, in order: (1) static text and (2) command output.
* Command execution when the panel widget is clicked, with separate commands for left-, middle-, and right-click actions. 

~~~
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 3: my Executor extension settings.</b></caption>
    <tr>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Field</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Left widget</b>
        </td>
        <td style="font-size: 20px; padding: 10px; text-align: center;">
            <b>Centre widget</b>
        </td>
        <td style="font-size: 20px; padding: 10px;">
            <b>Right widget</b>
        </td>
    </tr>
    <tr>
    <td style="font-size: 16px; padding: 10px;">
    <b>Index in widget</b>
    </td>
    <td style="font-size: 16px; padding: 10px;">3
    </td>
    <td style="font-size: 16px; padding: 10px;">2
    </td>
    <td style="font-size: 16px; padding: 10px;">2
    </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Output command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/left_widget_command.sh" target="_blank"><code>~/lfs_scripts/left_widget_command.sh</code></a> &mdash; displays the boot time and age of the system. The age is displayed as days/minutes/years hours:minutes:seconds. In my set up, it is used to generate output for the left widget. Runs every 60ms.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/centre_widget_command.sh" target="_blank"><code>~/lfs_scripts/centre_widget_command.sh</code></a> &mdash; displays CPU, RAM and root filesystem usage percentage and the number of the currently shown wallpaper / the total number of wallpapers in <code>~/wallpapers</code>. Runs every millisecond.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/updates_no.sh" target="_blank"><code>~/lfs_scripts/updates_no.sh</code></a> &mdash; checks for updates using the <a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/lfs-updates.sh" target="_blank"><code>updates</code></a> command in the shell profile. It displays <code>$in_progress󰔚 $updates_med_duration ($updates_iqr_duration)  $mod_time  $no_updates 󰂕 $no_missing_total  $no_failed$failed_version</code> where <code>$in_progress</code> is replaced with nothing if the <code>updates</code> command is not running, and <code>󰦕 ${percent}% </code> otherwise, where <code>$percent</code> is an approximation of how far through the running of <code>updates</code> we are. <code>$updates_med_duration</code> is the median duration, in minutes and seconds, of the run of <code>updates</code> based on <code>~/logs/updates_duration.log</code>. <code>$updates_iqr_duration</code> is the interquartile range of the duration of <code>updates</code> runs. <code>$mod_time</code> is replaced with the time the <code>updates</code> command last stopped running. <code>$no_updates</code> is replaced with the number of available package updates. <code>$no_missing_total</code> is replaced with the number of packages with missing inventories. <code>$no_failed</code> is replaced with the number of package versioning failures. <code>$failed_version</code>, if <code>~/log/failed_versioning.log</code> is not empty, is replaced by F and the number of packages with version failures in <code>~/logs/failed_versioning.log</code>. <code>updates</code> runs every 5 minutes - the median duration of <code>updates</code> runs - the IQR of duration of <code>updates</code> runs. <code>updates_no.sh</code> is run every millisecond.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Left-click command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-previous.sh" target="_blank"><code>~/lfs_apps/cycle-wallpaper-previous.sh</code></a> &mdash; show previous wallpaper.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic <a href="https://github.com/fusion809/lfs_scripts/blob/master/list-wallpapers.sh" target="_blank">~/lfs_scripts/list-wallpapers.sh</a></code> &mdash; displays the list of wallpapers in `~/wallpapers` with the currently shown wallpaper highlighted and centred.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "<a href="https://github.com/fusion809/NixOS-configs/blob/26.05/shell/user/21-lfs.sh" target="_blank">updatec</a>; exec zsh"</code> &mdash; updates the system's packages, including those installed via book instructions, custom packages and pip-managed packages and removes unneeded files. 
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Middle-click command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper-shuffle.sh" target="_blank"><code>~/lfs_apps/cycle-wallpaper-shuffle.sh</code></a> &mdash; show a randomly-selected wallpaper.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-extensions prefs executor@raujonas.github.io</code> &mdash; opens the settings dialog for Executor.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "source <a href="https://github.com/fusion809/lfs_scripts/blob/master/updates_no_func.sh" target="_blank">~/lfs_scripts/updates_no_func.sh</a>; silent_updates"</code> &mdash; runs <code>updates</code> to update the output shown in the widget.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Right-click command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_apps/blob/master/cycle-wallpaper.sh" target="_blank"><code>~/lfs_apps/cycle-wallpaper.sh</code></a> &mdash; show next wallpaper.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/open-wallpaper.sh" target="_blank"><code>~/lfs_scripts/open-wallpaper.sh</code></a> &mdash; opens the displayed wallpaper in Eye of GNOME.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <code>gnome-terminal -- zsh -ic "tail -f ~/updates.log"</code> &mdash; opens a terminal and follows the output of the <code>updates</code> command being used to generate the widget content.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Tooltip text</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Left click: previous wallpaper (Win+Z).<br/>
Middle click: shuffle wallpaper (Win+S).<br/>
Right click: next wallpaper (Win+W).<br/>
Win+N: show wallpaper whose number you will be asked to specify.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Left click: list wallpapers with displayed wallpaper centred and highlighted.<br/>
Middle click: open Executor settings (Win+E).<br/>
Right click: open wallpaper in EOG.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            Left click: run `update`.<br/>
Middle click: update notifications.<br/>
Right click: show log of last update check.
        </td>
    </tr>
    <tr>
        <td style="font-size: 16px; padding: 10px;">
            <b>Tooltip command</b>
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/left_widget_tooltip_command.sh" target="_blank"><code>~/lfs_scripts/left_widget_tooltip_command.sh</code></a> &mdash; generates a line describing the version of LFS/BLFS installed, along with the number of packages installed via different means, and package inventory commit numbers in a similar format as in the Fastfetch output. Also includes lines indicating how far into the current run of <code>autobuild &lt;package&gt;</code> the system is. 
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/centre_widget_tooltip_command_wrap.sh" target="_blank"><code>~/lfs_scripts/centre_widget_tooltip_command_wrap.sh</code></a> &mdash; lists selected wallpaper (indicated with <code>></code>) and the 25 wallpapers before and after this one. If there are not 25 wallpapers before the current one, it will show some of the last wallpapers in the collection before the wallpaper numbered 1 to ensure that 51 wallpapers are listed (including the one set as the desktop background). If there are not 25 wallpapers after the current one, it will show some of the first wallpapers in the collection after the final one in the list to ensure that 51 wallpapers are listed in total.
        </td>
        <td style="font-size: 16px; padding: 10px;">
            <a href="https://github.com/fusion809/lfs_scripts/blob/master/update-table.sh" target="_blank"><code>~/lfs_scripts/update-table.sh</code></a> &mdash; generates a more compact table of packages with updates, missing inventories and versioning failures.
        </td>
    </tr>
</table>
~~~

~~~
<br/>
<table style="border-collapse: collapse;">
    <caption style="font-size: 24px; padding: 10px; text-align: left;"><b>Table 4: Example tooltip contents.</b></caption>
    <tr>
    <td style="font-size: 20px; padding: 10px; text-align: center;">
    <b>Left</b>
    </td>
    <td style="font-size: 20px; padding: 10px; text-align: center;">
    <b>Centre</b>
    </td>
    <td style="font-size: 20px; padding: 10px; text-align: center;">
    <b>Right</b>
    </td>
    </tr>
    <tr>
    <td style="font-size: 16px; padding: 10px;">
    <img src="https://fusion809.github.io/images/executor-raujonas.github.io/Left_tooltip.png"/>
    </td>
    <td style="font-size: 16px; padding: 10px;">
    <img src="https://fusion809.github.io/images/executor-raujonas.github.io/Centre_tooltip.png"/>
    </td>
    <td style="font-size: 16px; padding: 10px;">
    <img src="https://fusion809.github.io/images/executor-raujonas.github.io/Right_tooltip.png"/>
    </td>
    </tr>
</table>
~~~

## Installing extensions via my web browser
BLFS did not provide a `gnome-browser-connector` package, which is required for installing GNOME extensions within one's browser. Manually compiling and installing it was fairly easy, however. That being said, whenever I tried to install an extension using it, I noticed that the extension was not successfully installed despite there being a folder in `~/.local/share/gnome-shell/extensions` for it. As this folder would be completely empty. Why? Well, running strace on GNOME shell revealed the problem was actually that the `gnome-browser-connector` was running `unzip` commands that assumed that Info-ZIP's unzip command was installed, not the `bsdunzip` variety provided by libarchive (which is the only one provided by BLFS or LFS). 

I tried compiling Info-ZIP's unzip, such as by following some [old BLFS instructions](https://www.linuxfromscratch.org/blfs/view/cvs/general/unzip.html) but this failed as Info-ZIP's unzip has not been updated since ~2009 and requires multiple intricate patches to get it to compile. The consolidated patch provided by BLFS was not even sufficient, even after I located the patch (the link provided in the book entry shared is actually dead, so I had to find a link to the patch elsewhere by Googling). 

Luckily, ChatGPT provided a script version of `unzip` that would run `bsdtar` in the background and could take all the arguments that `gnome-browser-connector` provided it. I have since included this script in my [custom package for libarchive](https://github.com/fusion809/lfs_packaging/tree/master/libarchive).
