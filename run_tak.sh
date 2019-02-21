#!/bin/sh

export PATH=/usr/local/bin:$PATH
export PATH=/usr/local/jdk1.8/bin:$PATH
export PATH=/lab/solexa_weng/testtube/trinityrnaseq-Trinity-v2.8.4:$PATH
export PATH=/lab/solexa_weng/testtube/TransDecoder-TransDecoder-v5.5.0:$PATH
export PATH=/lab/solexa_weng/testtube/transrate-1.0.3-linux-x86_64:$PATH
export PATH=/lab/solexa_weng/testtube/signalp-4.1:$PATH
export PATH=/lab/solexa_weng/testtube/targetp-1.1:$PATH
export PATH=/lab/solexa_weng/testtube/tmhmm-2.0c/bin:$PATH
export PATH=/lab/solexa_weng/testtube/ncbi-blast-2.8.1+/bin:$PATH
export PATH=/lab/solexa_weng/testtube/interproscan:$PATH
export PATH=/lab/solexa_weng/testtube/miniconda3/bin/:$PATH

./run_LSF.sh "$@"

