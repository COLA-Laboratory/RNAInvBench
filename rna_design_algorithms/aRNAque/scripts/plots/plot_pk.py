import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb


def get_bp_position(structure) :
    pk_pairs,cbrk, bp, nbp= [] , [], [], []

    pairs = {
        'bp' : [],
        'pk' : [],
        'nbp' : []
    }
    for i,elt in enumerate(structure) :
        if elt =='[' :
            pk_pairs.append(i)
        elif elt ==']' :
            pairs['pk'] += [(pk_pairs.pop(),i)]
        elif elt == '(' :
            bp += [i]
        elif elt == ')' :
            pairs['bp'] += [(bp.pop(),i)]
        elif elt == '{' :
            cbrk += [i]
        elif elt == '}' :
            pairs['pk'] += [(cbrk.pop(),i)]
        else :
            pairs['nbp'] += [i]
    return pairs

def bp_distance(target, strc2) :
    pos1 = get_bp_position(target)
    pair1 = set(pos1['bp'] + pos1['pk'])
    pos2 = get_bp_position(strc2)
    pair2 = set(pos2['bp'] + pos2['pk'])

    return len(pair1) + len(pair2) - 2*len(pair1.intersection(pair2))



def main() :
    #Load pk data
    ipknot_df = pd.read_csv("../../data/PseudoBase++/result_ipknot.csv")
    new_data = []
    s = 20
    labels= { str(t)+ " - " +str(t+s) : (t, t+s, []) for t in range(min(ipknot_df['Length'].values),max(ipknot_df['Length'].values),s)}

    print(labels)

    for key in labels.keys() :
        for l in range(labels[key][0],labels[key][1]):
            val = ipknot_df[ipknot_df["Length"]==l].values.tolist()
            if len(val) > 0 :
                for elt in val :
                    nl = list(elt)
                    nl [-1] = nl[-1] + 1
                    new_data += [nl+[key]+[bp_distance(nl[2], nl[3])]]
    print(len(ipknot_df), len(new_data))
    new_dt = pd.DataFrame(new_data, columns=["id","sequence","structure","target", "Length","Hamming distance",'bp_density'
    ,'Mutation mode','Number of generations', "Length group", "BP distance"])
    op_df = new_dt[new_dt["Mutation mode"]=='OP']
    op_med = {}



    levy_df = new_dt[new_dt["Mutation mode"]=='Levy']
    levy_med = {}
    for key in labels.keys():
        gens = op_df[op_df["Length group"]==key]["Number of generations"].values
        print("One point: "+key, np.median(gens))

        op_med[key] = np.median(gens)
        gens = levy_df[levy_df["Length group"]==key]["Number of generations"].values
        print("Levy: "+key, np.median(gens))
        levy_med[key] = np.median(gens)


    figure = plt.figure(constrained_layout=True, figsize=(9,6))
    gs = figure.add_gridspec(nrows=2, ncols=1, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    #sb.set(rc={'figure.figsize': (20., 8.27)})
    #plt.ylabel('Hamming Distance')
    plt.title("(A)")
    plt.xticks(rotation=90, fontsize=12)
    plt.title("(A)", fontsize=15)
    ax.set_ylabel('BP distance',fontsize=12)
    #ax.set_xlabel('Length Group', fontsize=12)
    sb_bx = sb.boxplot(ax=ax, y='BP distance', x='Length group', hue='Mutation mode', data=new_dt, palette={
    "Levy": "deepskyblue",
    "OP": "darkorange"
    })
    sb_bx.set(xticklabels=[])
    sb_bx.set(xlabel=None)
    handles, _ = sb_bx.get_legend_handles_labels()          # Get the artists.
    sb_bx.legend(handles, ["Lévy mutation", "Local mutation"], loc="best") # Associate manually the artists to a label.


    ax = figure.add_subplot(gs[1,0])
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #plt.xticks(rotation=90)
    plt.title("(B)", fontsize=15)
    ax.set_ylabel('BP distance', fontsize=12)
    ax.set_xlabel('Length group', fontsize=12)
    ax.set(yscale="log")
    sb_ax = sb.boxplot(ax=ax, y='Number of generations', x='Length group', hue='Mutation mode', data=new_dt,palette={
    "Levy": "deepskyblue",
    "OP": "darkorange"
    })

    plt.legend([],[], frameon=False)
    plt.savefig('../../images/PseudoBase++/fig4.pdf')
    plt.show()

    diff = [op_med[k]-levy_med[k] for k in op_med.keys()]
    print(diff, np.mean(diff))

    diff = [op_med[k]-levy_med[k] for k in op_med.keys()]
    print(diff, np.median(diff))
    figure = plt.figure(constrained_layout=True, figsize=(7,5))
    gs = figure.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.48, wspace=0.05)
    ax = figure.add_subplot(gs[0,0])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(A)",fontsize=15)
    plt.xlabel(r"Median number of generations ($t$)")
    plt.hist(levy_med, bins=10, label="Levy", color="deepskyblue")
    plt.hist(op_med, bins=10, alpha=0.3,label="OP", color="darkorange")
    #plt.legend()

    op_success = [(20-list(gens).count(200))/0.2 for gens in op_df]
    levy_success = [(20-list(gens).count(200))/0.2 for gens in levy_df]
    ax = figure.add_subplot(gs[0,1])

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.title("(B)",fontsize=15)
    plt.xlabel("Success rate (%)")
    plt.hist(levy_success, bins=10, label="Levy", color="deepskyblue")
    plt.hist(op_success, bins=10, alpha=0.4,label="OP", color="darkorange")
    plt.legend()
    plt.savefig("../../images/PseudoBase++/pk_histo2.pdf")
    plt.show()

    print()



if __name__=="__main__" :
    main()
