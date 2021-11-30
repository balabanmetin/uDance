
include: "workflows/decompose.smk"

rule all:
    input: "test.txt"

rule clean:
    shell:
        """
            rm -r test.txt
        """
