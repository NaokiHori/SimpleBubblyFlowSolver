def add(macros):
    # number of grid points
    macros["nx"] = "{N_{\\vx}}"
    macros["ny"] = "{N_{\\vy}}"
    macros["nz"] = "{N_{\\vz}}"
    # summation symbols for different locations
    macros["sumxf"] = "\\sum_{i = \\frac{1}{2}}^{\\nx + \\frac{1}{2}}"
    macros["sumxc"] = "\\sum_{i = 1}^{\\nx}"
    macros["sumyf"] = "\\sum_{j = \\frac{1}{2}}^{\\ny - \\frac{1}{2}}"
    macros["sumyc"] = "\\sum_{j = 1}^{\\ny}"
    macros["sumzf"] = "\\sum_{k = \\frac{1}{2}}^{\\nz - \\frac{1}{2}}"
    macros["sumzc"] = "\\sum_{k = 1}^{\\nz}"
    # discrete momentum balance
    macros["dmomadv"] = [
            "-"
            "\\frac{1}{J}"
            "\\dif{"
            "  \\left("
            "    \\ave{"
            "      \\frac{J}{h_{\\gcs^{#2}}}"
            "      \\rho u_{#2}"
            "    }{\\gcs^{#1}}"
            "    \\ave{u_{#1}}{\\gcs^{#2}}"
            "  \\right)"
            "}{\\gcs^{#2}}"
            , 2
    ]
    macros["dmompre"] = [
            "-"
            "\\frac{1}{h_{#1}}"
            "\\dif{p}{\\gcs_{#1}}"
            , 1
    ]
    macros["dmomdif"] = [
            "+"
            "\\frac{1}{J}"
            "\\dif{}{\\gcs^{#2}}"
            "\\left("
            "   \\frac{J}{h_{\\gcs^{#2}}}"
            "   \\tau_{#1 #2}"
            "\\right)"
            , 2
    ]

