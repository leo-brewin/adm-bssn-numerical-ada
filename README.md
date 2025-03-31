# Evolving initial data using the ADM and BSSN equations. II.

The companion [repo][1] provides a series of preprocessing tols and a set of Ada templates to build a complete Ada code (to model the evolution of a universe). The starting point is a set of mathematical equations which are translated (via third party tools) into fragments of C-code. These fragments are then merged with a set of Ada templates to produce the final Ada code. The journey is non-trivial and it does require a number of 3rd party tools to be installed ([Cadabra][2] to read and process the mathematics and [sympy][3] for the code generation). That may be a bridge too far for those interested in just the final Ada code. And that is the purpose of the repo -- to serve the final Ada code without the need for any third party tools.

## Running the codes

It's really very simple, just use `make` to start the ball rolling.

```sh
$ make
```

If all goes well (as it should), you will see two short summaries of results on the terminal. These summaries are also saved as `adm/code/adm-history.txt` and `bssn/code/bssn-history.txt` which you can compare against the expected results in the `exepcted/` directory. There should also be a new `adm-bssn-plots.pdf` file (which should match exactly its counterpart in the `exepected/` directory).

There are a number of other targets in the Makefile that you might find useful (e.g., `make code` to just compile but not run the codes).

In the unlikely event that the read/write/execute premissions of one or more of the scripts is incorrect, you can set the correct permissions using

```sh
$ source CHMOD.txt
```

## Tinkering

There are four Ada programs, two to set some initial data (`adm/code/adminitial.adb` and `bssn/code/bssninitial.adb`) and two to evolove that data (`adm/code/admevolve.adb` and `bssn/code/bssnevolve.adb`).

Each program can be compiled and run using a like named script. So to compile and run the `admevolve.adb` code, you can use

```sh
$ (cd adm/code; admevolve.sh)
```

Each program accepts a number of command line options (all long form options). These can be passed directly to the binary or as arguments to the matching shell script.

For a list of the options and their expected parameters, you can use, for example,

```sh
$ (cd adm/code; admevolve.sh --Help)
```

## License

All files in this collection are distributed under the [MIT][10] license. See the file LICENSE.txt for the full details.

 [1]: https://github.com/leo-brewin/adm-bssn-numerical
 [2]: https://cadabra.science
 [3]: https://www.sympy.org/
[10]: https://opensource.org/licenses/MIT
