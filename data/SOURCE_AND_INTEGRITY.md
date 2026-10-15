# CHC data: source and integrity

The application uses the estimation dataset distributed with:

Martha J. Bailey and Andrew Goodman-Bacon (2015), "The War on Poverty's
Experiment in Public Medicine: Community Health Centers and the Mortality of
Older Americans," *American Economic Review* 105(3), 1067–1104.
[Article DOI](https://doi.org/10.1257/aer.20120070).

The original replication archive was available from
`https://www.aeaweb.org/aer/data/10503/20120070_data.zip` and is catalogued as
openICPSR project E112871V1,
[replication archive DOI](https://doi.org/10.3886/E112871V1).
The working project obtained an unchanged copy of that AEA ZIP through the
Internet Archive. This package includes its unchanged estimation file and
original documentation, rather than a second copy of the ZIP.

| Included file | Description |
| --- | --- |
| `aer_data.dta` | Original estimation data; read without modification |
| `original_readme.pdf` | Documentation supplied with the original archive |

SHA-256 of `aer_data.dta`:

```text
370e38b04009750dc444d8a1f9a9bcd1105c1351174ee5c08772e378e35b3b1d
```

For source tracking, the original ZIP's SHA-256 is
`ce956fce08bce77ea463725a52cafe63a8a77f5f1b7ea624a7d4b6e9563d53c0`.
The full application verifies the DTA checksum before reading it. No network
access, additional data purchase, or manual sample preparation is needed.
Sample restrictions and variable transformations are in `code/application.R`;
`results/application/sample_audit.csv` records the sample construction.

These materials retain the original authors' rights and applicable source
terms. No new license for the data or original documentation is asserted here.
