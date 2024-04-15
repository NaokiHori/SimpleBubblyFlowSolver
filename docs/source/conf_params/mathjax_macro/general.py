def add(macros):
    # coordinate
    macros["vx"] = "{x}"
    macros["vy"] = "{y}"
    macros["vz"] = "{z}"
    # a symbol used to denote general coordinate
    macros["gcs"] = "{\\xi}"
    macros["gx"] = "{\\gcs^{\\vx}}"
    macros["gy"] = "{\\gcs^{\\vy}}"
    macros["gz"] = "{\\gcs^{\\vz}}"
    # velocity
    macros["ux"] = "{u_{\\vx}}"
    macros["uy"] = "{u_{\\vy}}"
    macros["uz"] = "{u_{\\vz}}"
    # quadratic quantities
    macros["kx"] = "{k_{\\vx}}"
    macros["ky"] = "{k_{\\vy}}"
    macros["kz"] = "{k_{\\vz}}"
    # scale factors
    macros["hx"] = "{h_{\\gx}}"
    macros["hy"] = "{h_{\\gy}}"
    macros["hz"] = "{h_{\\gz}}"
    # jacobian determinant divided by the scale factors
    macros["jhx"] = "{\\frac{J}{\\hx}}"
    macros["jhy"] = "{\\frac{J}{\\hy}}"
    macros["jhz"] = "{\\frac{J}{\\hz}}"
    # differentiations
    macros["pder"] = ["{\\frac{\\partial #1}{\\partial #2}}", 2]
    macros["dder"] = ["{\\frac{\\delta #1}{\\delta #2}}", 2]
    # discrete operators
    macros["ave"] = ["{\\overline{#1}^{#2}}", 2]
    macros["dif"] = ["{\\delta_{#2} {#1}}", 2]
    macros["vat"] = ["{\\left. {#1} \\right|_{#2}}", 2]
