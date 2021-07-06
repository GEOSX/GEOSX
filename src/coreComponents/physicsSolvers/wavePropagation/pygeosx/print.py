def print_group(group, indent=0):
    print("{}{}".format(" " * indent, group))

    indent += 4
    print("{}wrappers:".format(" " * indent))

    for wrapper in group.wrappers():
        print("{}{}".format(" " * (indent + 4), wrapper))
        print_with_indent(str(wrapper.value(False)), indent + 8)

    print("{}groups:".format(" " * indent))

    for subgroup in group.groups():
        print_group(subgroup, indent + 4)
        
        
def print_with_indent(msg, indent):
    indent_str = " " * indent
    print(indent_str + msg.replace("\n", "\n" + indent_str))
    

def print_flag(shot_list):
    i = 0
    for shot in shot_list:
        print("Shot " + str(i) + " status : " + shot.getFlag())
        i += 1
    print("\n")
     
def print_pressure(pressure, ishot):
    print("\n" + "Pressure value at receivers for configuration " + str(ishot) + " : \n")
    print(pressure)
    print("\n")
    
def print_shot_config(shot_list, ishot):
    print("\n \n" + "Shot configuration number " + str(ishot) + " : \n") 
    print(shot_list[ishot])
