# `~/logs`
`~/logs` contains assorted log files for the LFS system. Although, it is important to note that the logs from building packages are stored in `~/build_logs` and due to how large these logs are, they are not git controlled. 

* `backup_updates_duration.log` &mdash; a backup of an old version of `updates_duration.log`. 
* `book_hash.log` &mdash; 7-character git hash for the `/var/lib/book-packages` git repo. 
* `custom_hash.log` &mdash; 7-character git hash for the `/var/lib/custom-packages` git repo. 
* `firefox.log`, `gh-bin.log`, `julia-bin.log` and `rust-bin.log` &mdash; these document the version and size (in human-readable format) of these packages when they were last installed. These are used to produce `pkgs_table` output. 
* `inventory_commit_no_long.log` &mdash; contains `Package inventory commit number:    585`, for instance, where 585 is the total number of `/var/lib/custom-packages` git commits.
* `lfs-time.log` &mdash; time that the `lfs-downloader.sh` script last ran.
* `os_version.log` &mdash; version string of LFS virtual machine `󰌽 r13.1-17,284`.
* `packages_no.log` &mdash; package counts for virtual machine `852 [ 735,  2,  86,  29]`.
* `packages_no_long.log` &mdash; longer package counts string `852 [  735 (󰊢 585)  2  86  29]`. 
* `pkgs_by_alpha.log`, `pkgs_by_bd.log` and `pkgs_by_size.log` &mdash; contains a table displaying packages with their installed sizes, build times, versions and descriptions &mdash; sorted alphabetically by name, build time and size, respectively.
* `updates_duration.log` &mdash; contains how long `updates` runs used to update GNOME top panel bar have taken. 
* `updates.log` &mdash; the output of the most recent run of `updates`.