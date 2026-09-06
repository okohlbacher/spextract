# Security policy

## Reporting a vulnerability

Please report security issues privately via GitHub's "Report a vulnerability" button on the
Security tab, rather than opening a public issue.

## Scope

SpeXtractor is a batch command-line tool. It reads instrument data and writes mzML; it opens no
sockets, runs no server, and requires no credentials. The realistic security surface is therefore
**untrusted input files**: a malformed or hostile Bruker `.d`, `.mzML` or `.mzpeak` reaching the
parser. Reports of crashes, out-of-bounds reads or unbounded allocation triggered by an input file
are in scope and welcome.

Most parsing is done by OpenMS and by the vendor libraries. If the flaw is there, we will help
route it upstream.

## Data handling

SpeXtractor reads local files and writes local files. It does not phone home, collect telemetry, or
transmit anything. If you find otherwise, that is a security bug — report it.
