# Type IIS Digest & Ligation Checker

A static GitHub Pages application for checking Type IIS restriction digests and ligation compatibility.

## What it does

- Accepts FASTA, GenBank, or plain DNA sequence input
- Scans both strand orientations for selected Type IIS restriction sites
- Displays each hit as double-stranded DNA with cut positions shown explicitly
- Supports a restriction-site block builder for typed sequence design
- Enumerates simple sticky or blunt ends for ligation compatibility testing
- Runs entirely in the browser, so GitHub Pages is enough

## Files

- `index.html` — app shell and UI
- `styles.css` — layout and styling
- `enzymes.js` — curated Type IIS enzyme catalog
- `app.js` — sequence parsing, scanning, cut mapping, and compatibility logic

## Deploy to GitHub Pages

1. Create a new GitHub repository.
2. Upload these files to the repository root.
3. Commit to the default branch.
4. In GitHub, go to **Settings → Pages**.
5. Under **Build and deployment**, choose **Deploy from a branch**.
6. Select your default branch and the `/ (root)` folder.
7. Save. GitHub Pages will publish the app.

## Notes

- The enzyme catalog is curated from the current NEB Type IIS selection chart.
- Simple cutters are fully supported for ligation compatibility checks.
- Dual-flank cutters are detected and their cuts are displayed, but they are not enumerated into the automatic ligation pair scanner in this first version.
- If you want full support for all REBASE Type IIS/IIC edge cases, the catalog and digest model can be extended further.
