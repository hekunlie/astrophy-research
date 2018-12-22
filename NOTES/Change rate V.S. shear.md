Change rate V.S. shear

Profile: Disc-dominated & Bulge-dominated

Ellipticity: 0, 0.3, 0.5, 0.8

BTR: 0, 0.2, 0.4, 0.6, 0.8, 1.0

Radius (arcsec): 0.2, 0.4, 0.6, 0.8, 1.0

SEX filter: gauss2.0, gauss3.0, gauss4.0

Threshold: 1.5$\sigma$, 2$\sigma$ 





```mermaid
graph TD
Profile --> B[Disc-dominated]
B[Disc-dominated] --> B1[Elliptical]
B1 --> B11[...]
B1 --> B12[BTR_i]

B12 --> B121[....]
B12 --> B122[Radius_i]

B --> B2[Circle]
B2 --> B21[...]


Profile --> C[Bulge-dominated]
C --> C1[...]

```

