## Annotation of the peptide search database
Generate file list and their MD5 checksums under `intermediates` by:

```bash
fd --type file . intermediates \
    | sort \
    | parallel --keep-order -j6 --bar 'gmd5sum {}' \
    > intermediates_files_md5_checksum.txt
```