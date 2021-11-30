
rule decompose:
    output: "test.txt"
    shell:
        """
            touch test.txt
        """
