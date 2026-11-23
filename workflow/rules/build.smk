# build stamp + frozen-input verification (first rules of every build)
rule stamp:
    output: STAMP
    log: LOG("stamp")
    shell: WC("stamp") + " > {log} 2>&1"


rule frozen:
    input: STAMP
    output: FROZEN_OK
    log: LOG("frozen")
    shell: WC("frozen") + " > {log} 2>&1"
