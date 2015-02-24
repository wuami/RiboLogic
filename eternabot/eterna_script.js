// -------------------------------------------------------------
// Things that should go in the Library class (imo)
// -------------------------------------------------------------

var sid = "3386938";
var dbgwin = window.open("","EternaScript_"+sid,"width=400,height=800,toolbar=no,location=no,directories=no,status=no,menubar=no,scrollbars=yes,resizable=yes");
var dbgdoc = dbgwin.document;

dbgdoc.title = "EternaScript ID:"+sid+" output";

var dbg_out = function(txt, elemtype) {
    var elem = dbgdoc.createElement(elemtype);
    var x = dbgdoc.createTextNode(txt);
      
    elem.appendChild(x);
    dbgdoc.body.appendChild(elem);
};

// -------------------------------------------------------------

var status = "";

var applet = document.getElementById('maingame');
// console.log(applet);
if (applet == null) {
    status = "maingame missing";
    return status;
}

var seq = applet.get_sequence_string();
// console.log(seq);
dbg_out(seq, "div");

var locks = applet.get_locks();
// console.log(locks);
dbg_out(locks.toString(), "div");

var targets = applet.get_targets();
// console.log(targets);
dbg_out(JSON.stringify(targets, undefined, 2), "pre");

var constraints = applet.get_constraints();
// console.log(constraints);
dbg_out(JSON.stringify(constraints, undefined, 2), "pre");

applet.set_script_status("wait for input");
var niter = prompt("How many iterations would you like to run?", 1000)
if (niter === null) return "cancelled";
niter = parseInt(niter,10)
if (niter !== NaN or niter <= 0) return "invalid input";

var index_array = new Array(); 
for ( i = 0; i < locks.length; i++ ) {
    if (!locks[i]) {
        index_array.push(i);
    }
}

function random_base() {
    return Lib.bases[Math.floor(Math.random() * 4)];
}

function mutate(sequence) {
    rindex = index_array[Math.floor(Math.random() * index_array.length)];
    var seq_array = seq.split("");
    seq_array[rindex] = random_base();
    return seq_array.join("");
}

function bp_distance(secstruct1, secstruct2, constraints) {
    if (typeof constraints === "undefined") {
        constraints = Array.apply(null, new Array(secstruct1.length)).map(Number.prototype.valueOf,0);
    }

    pairmap1 = getPairmap(secstruct1);
    pairmap2 = getPairmap(secstruct2);

    var dist = 0;
    for (i = 0; i < constraints.length; i++ ) {
        if (constraints[i]) {
            if (pairmap1[i] != pairmap2[i]) {
                if (pairmap1[i] > i) {
                    dist++;
                }
                if (pairmap2[i] > i) {
                    dist++;
                }
            }
        }
    }
    return dist;
}

function score_sequence(sequence) {
    score = 0;
    for (target in targets) {
        if (target['type'] == "single") {
            score += bp_distance(target['secstruct'],Lib.fold(sequence));
        } else if (target['type'] == "oligo") {
            secstruct = Lib.fold(target['oligo_sequence'] + sequence);
            score += bp_distance(target['secstruct'],secstruct,target['structure_constraints']);
            if ('anti_secstruct' in target) {
                score += bp_distance(Array(target['anti_secstruct'].length).join("."),secstruct,target['anti_structure_constraints']);
            }
        }
    }
    return score;
}

current_seq = seq;
current_score = score_sequence(current_seq);
//best_seq = current_seq;
//best_score = current_score;
i = 0;
for ( ; i < niter; i++ ) {
    mut_seq = mutate(current_seq);
    mut_score = score_sequence(mut_seq);
    if (mut_score <= current_score) {
        current_seq = mut_seq;
        current_score = current_seq;
        var ok = applet.set_sequence_string(new_seq);
        dbg_out("mutation " + (ok? "succeeded" : "FAILED")  , "div");
        if (current_score == 0) {
            break;
        }
    }
}


status = "done: " + i + " iterations";
return status;
