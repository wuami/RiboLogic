import os
import sys
import subprocess
import re

class Varna(object):
    def __init__(self):
        self.varna_command = "java -cp ./resources/VARNA.jar fr.orsay.lri.varna.applications.VARNAcmd"


    def new_image(self,motif,filename,colormapstyle="red",highlight_res=[],highight_bp=[]):
        sstructure,sequence,design_sequence = motif.generate_secondary_structure()

        #remove last + sign
        sstructure_str = "".join(sstructure)[:-1]
        sequence_str = "".join(sequence)[:-1]

        colormap = [0 for x in range(len(sequence_str))]

        for r in highlight_res:
            index = self.get_res_pos(motif,r)
            colormap[index] = 1


        colormap_str = ";".join([str(x) for x in colormap])

        subprocess.call(self.varna_command + " -sequenceDBN " + sequence_str + " -structureDBN \"" + sstructure_str + "\" -o " + filename + " -colorMap \"" + colormap_str + "\" -colorMapStyle " + colormapstyle, 
            shell=True)

    def new_image_by_str(self,filename,sstructure_str,sequence_str,colormapstyle="red",highlight_res_num=[],size=[],colormap_str=""):
        

        #colormap = [0 for x in range(len(sequence_str))]

        #for i in highlight_res_num:
        #   #print len(colormap),i
        #   pos = i-1 
        #   nplus = 0

        #   for i,e in enumerate(sequence_str):
        #       if i == pos+nplus:
        #           break

        #       if e == "+":
        #           nplus += 1

        #   colormap[pos+nplus] = 1

        #colormap_str = ";".join([str(x) for x in colormap])

        #subprocess.call(self.varna_command + " -sequenceDBN " + sequence_str + " -structureDBN \"" + sstructure_str + "\" -o " + filename + " -colorMap \"" + colormap_str + "\" -colorMapStyle " + colormapstyle +  " -bp \"#000000\" -periodNum 10 -spaceBetweenBases 0.7", 
        #   shell=True)

        subprocess.call(self.varna_command + " -sequenceDBN " + sequence_str + " -structureDBN \"" + sstructure_str + "\" -o " + filename + " -colorMap \"" + colormap_str + "\" -flat false", 
            shell=True)

        return

        f = open(filename)
        lines = f.readlines()
        f.close()

        p = re.compile(r"=\"(\d+\.\d+)")
        max_x = 0
        max_y = 0

        for l in lines:
            spl = re.split("\s+",l)
            if len(spl) < 5:
                continue

            if spl[0] != "<line":
                continue

            for i in range(1,5):
                if p.search(spl[i]):
                    m = p.search(spl[i])
                    num = float(m.group(1))
                    if i == 1 or i == 3:
                        if num > max_x:
                            max_x = num
                    else:
                        if num > max_y:
                            max_y = num
            
        lines[4:5] = "<svg width=\"300pt\" height=\"300pt\" viewBox=\"0 0 "+str(max_x+20)+" "+str(max_y+20)+"\" xmlns=\"http://www.w3.org/2000/svg\">\n"

        del lines[-84:]

        f = open(filename,"w")
        for l in lines:
            f.write(l)

        f.write("</svg>")
        f.close()


    def get_res_pos(self,motif,res):

        for i,c in enumerate(motif.structure.chains):
            index = 0
            if res not in c.residues:
                continue

            for j,c1 in enumerate(motif.structure.chains):
                if j <= i:
                    continue
                index += len(c1.residues)+1

            return index + c.residues.index(res)
        
