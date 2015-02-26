// -------------------------------------------------------------
// Things that should probably go in the Library class (imo)
// -------------------------------------------------------------

var sid, dbgwin, dbgdoc;
var dbg = false;

if (dbg) {
  sid = "3387609";
  dbgwin = window.open("","EternaScript_"+sid,"width=400,height=800,toolbar=no,location=no,directories=no,status=no,menubar=no,scrollbars=yes,resizable=yes");
  dbgdoc = dbgwin.document;

  dbgdoc.title = "EternaScript ID:"+sid+" output";
}

var dbg_out = function(txt, elemtype) {
  if (!dbg) return;
  var elem = dbgdoc.createElement(elemtype);
  var x = dbgdoc.createTextNode(txt);
  
  elem.appendChild(x);
  dbgdoc.body.appendChild(elem);
};

// -------------------------------------------------------------

var applet = document.getElementById('maingame');
if (applet == null) return "maingame element missing";

var seq = applet.get_sequence_string();
// console.log(seq);
//dbg_out(seq, "div");

var locks = applet.get_locks();
// console.log(locks);
//dbg_out(locks.toString(), "div");

var targets = applet.get_targets();
// console.log(targets);
//dbg_out(JSON.stringify(targets, undefined, 2), "pre");

var constraints = applet.get_constraints();
// console.log(constraints);
//dbg_out(JSON.stringify(constraints, undefined, 2), "pre");

applet.set_script_status("wait for input");
var niter = prompt("How many iterations would you like to run?", 100);
if (niter === null) return "cancelled";
niter = parseInt(niter,10);
if (isNaN(niter) || niter <= 0) return "invalid input";

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
    var rindex = index_array[Math.floor(Math.random() * index_array.length)];
    var seq_array = sequence.split("");
    seq_array[rindex] = random_base();
  	return seq_array.join("");
}

function bp_distance(secstruct1, secstruct2, constraints) {
    if (typeof constraints === "undefined") {
        constraints = Array.apply(null, new Array(secstruct1.length)).map(Number.prototype.valueOf,1);
    }

  	var r1 = new RNA(secstruct1);
    var pairmap1 = r1.getPairmap(secstruct1);
    var r2 = new RNA(secstruct2);
    var pairmap2 = r2.getPairmap(secstruct2);
  
    var dist = 0;
    for (var ii = 0; ii < constraints.length; ii++ ) {
        if (constraints[ii]) {
          	if (pairmap1[ii] === undefined) {
            	pairmap1[ii] = -1;
          	}
          	if (pairmap2[ii] === undefined) {
            	pairmap2[ii] = -1;
          	}
            if (pairmap1[ii] != pairmap2[ii]) {
                if (pairmap1[ii] > ii) {
                    dist++;
                }
                if (pairmap2[ii] > ii) {
                    dist++;
                }
            }
        }
    }
    return dist;
}

function score_sequence(sequence) {
    var score = 0;
    for (var ii = 0; ii < targets.length; ii++) {
        target = targets[ii];
        if (target['type'] == "single") {
            var secstruct = applet.fold(sequence);
            score += bp_distance(target['secstruct'],secstruct);
        } else if (target['type'] == "oligo") {
            var secstruct = applet.fold(target['oligo_sequence'] + sequence);
            score += bp_distance(target['secstruct'],secstruct,target['structure_constraints']);
            if ('anti_secstruct' in target) {
                score += bp_distance(Array(target['anti_secstruct'].length).join("."),secstruct,target['anti_structure_constraints']);
            }
        }
    }
    return score;
}


function prob(dist, new_dist) {
    return Math.exp(-(new_dist-dist));
}

var current_seq = seq;
var current_dist = score_sequence(current_seq);
var best_seq = current_seq;
var best_dist = current_dist;
var mut_seq;
var mut_dist;
var i = 0;
for ( ; i < niter; i++ ) {
    mut_seq = mutate(current_seq);
    mut_dist = score_sequence(mut_seq);
    if (Math.random() <= prob(current_dist, mut_dist)) {
        current_seq = mut_seq;
        dbg_out(current_seq, "div");
        current_dist = mut_dist;
        dbg_out(current_dist, "div");
        if (current_dist < best_dist) {
            best_seq = current_seq;
            best_dist = current_dist;
            if (current_dist == 0) {
                break;
            }
        }
    }
}
dbg_out(best_seq, "div");

var ok = applet.set_sequence_string(best_seq);
dbg_out("mutation " + (ok? "succeeded" : "FAILED")  , "div");
 
var status = "done: " + i + " iterations";
return status;