# SingleNucSeq-Website


## Quicky Instalation Instructions

Check code out from github.  Go into the `cestaan` directory.

`mkdir data source_data scripts/source_data`, then copy data there.*
`cd scripts` and run the `create_database.py` command.  Move the
`scripts/data/output.sqlite` file to `data/output.sqlite`.

Check `app.py`, `generate_umap.R`, and `generate_violin.R` and make
sure `CESTAAN_ROOT` is logical.†

`mkdir venv`, and then `virtualenv venv`.  Use `venv/bin/pip` to
install pacages listed in `requirements.txt`.

To run it first do, in the root `cestaan` directory do
`venv/bin/activate`, and then `shiny run --port 8070 app.py`.  Or set
up a system service with a file like:

```
[Unit]
Description=Start CeSTAAN shiny thingy

[Service]
User=ctvm1
Group=ctvm1

WorkingDirectory=/var/www/cestaan
Environment=BASH_ENV=./venv/bin/activate
ExecStart=/usr/bin/bash -c "shiny run --port 8070 app.py"

[Install]
WantedBy=multi-user.target
```

The point your browser to port 8070, or set up a web server to proxy to that port.

### Footnotes

\* Not yet sure where to get the data.

† Should create global config.
