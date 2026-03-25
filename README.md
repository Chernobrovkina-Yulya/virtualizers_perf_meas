# Performance measurement
This code contain C linpack benchmark from  https://people.sc.fsu.edu/~jburkardt/c_src/linpack_bench/linpack_bench.html and
script to perform measurements on original and obfuscated programs.

Following virtualizers were tested:
- Softcom LLVM obfuscator
- Tigress
- VMProtect
- Themida

**Note:** Obfuscated binaries are not distributed with this repository.  
All third‑party obfuscators (Tigress, VMProtect, Code Virtualizer) are used under their respective demo licenses and must be downloaded separately from the official sources.

# Binary files preparation

Following functions were obfuscated during performance measurement:
```c
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
int dgefa ( double a[], int lda, int n, int ipvt[] );
void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job );
void dscal ( int n, double sa, double x[], int incx );
int idamax ( int n, double dx[], int incx );
double r8_random ( int iseed[4] );
double *r8mat_gen ( int lda, int n );
double r8mat_norm_li ( int lda, int m, int n, double A[] );
```

## Original binary
Compile source linpack file using clang:
```bash
clang -O0 linpack_bench.c -o linpack
```

## Softcom LLVM obfuscator
Use `llvm_linpack_bench.c` source file. It has comments for virtualization pass in all function above.

Compile C source file using clang version with obfuscation passes:
```bash
clang -O0 --target=x86_64_obf-unknown-linux-gnu llvm_linpack_bench.c -o llvm_linpack
```
**Note**: --target option is needed to activate virtualization pass.

## Tigress
Tigress is source to source C obfuscator, meaning result of obfuscation - another C source file.

Download Tigress deb package for Linux from official site https://tigress.wtf/ and install it.
```bash
 sudo dpkg --force-architecture -i file.deb
```
Put `tigress.h` header in current directory.

Add include of `tigress.h` header file in begging of linpack source code. Without header tigress won't obfuscate linpack.

Use following script to virtualize linpack using tigress:
```bash
tigress --Compiler=clang --Transform=Virtualize --Functions=comma_separated_functions_from_above \
        linpack_bench.c --out=tigress_linpack_bench.c
```
Compile result C file using clang:
```bash
clang -O0 tigress_linpack_bench.c -o tigress_linpack
```

## VMProtect (CLI version)
Install VMProtect Demo version for Linux from official site https://vmpsoft.com/ and install it.

functions `VMProtectBeginVirtualization` and `VMProtectEnd` define blocks of code to virtualize. `vmprotect_linpack_bench.c` source file already has this calls for functions above.

Compile vmprotect source file using following script:
```bash
clang -O0 -I/path_to/VMProtectSDK.h -L/path_to/libVMProtectSDK64.so vmprotect_linpack_bench.c -o vmprotect_linpack_tmp \
        -lVMProtectSDK64
```

Fed temporary binary file through vmprotect tool and get result `vmprotect_linpack` bin:
```bash
vmprotect_con vmprotect_linpack_tmp vmprotect_linpack
```

## Themida (Code Virtualizer)
**Note**: Themida and Code Virtualizer have the same virtual machines, so use Code Virtualizer for convenience.

Download Code Virtualizer Demo version from official site https://www.oreans.com/CodeVirtualizer.php.

**Note**: Code Virtualizer could only be run on Windows, but support ELF files obfuscation.

To use Code Virtualizer user should insert macros `VIRTUALIZER_START` and `VIRTUALIZER_END` in the beginning and end of each function to virtualize. `themida_linpack_benchmark.c` already has all macros and `VirtualizerSDK.h` header include.

Compile source file:
``` bash
clang -O0 -I/path_to_VirtulizerSDK.h themida_linpack_benchmark.c -o themida_linpack
```
Temporary ELF file would contain symbols for Code Virtualizer, extract file on Windows OS using any convenient method.

Launch Code Virtualizer on Windows and load `themida_linpack` ELF file. In Demo version only available option - use Falcon Tiny VM, so just click protect button. Result file by default will be called `themida_linpack_protected`, extract it back on Linux OS.

**NOTE**: for some reason result protected ELF would crash with SEgmentation Fault, so Themida results are not present in performance measurements.

# Measurements
execute python script:
```bash
python measurement_launcher.py
```

Enter matrix size and wait for results. csv file with result will be stored.

**Note**: for a 500 matrix size, VMProtect starts to slow down a lot.