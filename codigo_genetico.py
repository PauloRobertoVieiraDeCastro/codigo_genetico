class DNA:
    def __init__(self,sequencia):
        self._sequencia = sequencia.upper()

    @property
    def get_dna(self):
        return self._sequencia

    @staticmethod
    def codon_aminoacido():
        codon_tabela = {"UUU":"F","UUC":"F","UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L","AUC":"I","AUU":"I","AUA":"I","AUG":"M",
                 "GUU":"V","GUC":"V","GUA":"V","GUG":"V","UAC":"Y","UAU":"Y","UGU":"C","UGC":"C","UGG":"W","UCU":"S","UCC":"S","UCA":"S","UCG":"S"}
        return codon_tabela

    def __estabilidade(self):
        C = self._sequencia.count("C")
        G = self._sequencia.count("G")
        tamanho = len(self._sequencia)
        return 100*(C + G)/tamanho

    @property
    def get_estabilidade(self):
        r = self.__estabilidade()
        return r
    
class Transcricao(DNA):
    def __init__(self,sequencia):
        super().__init__(sequencia)    

    def transcriptase(self):
        rna = []
        for i in self._sequencia:
            if(i == "A"):
                rna.append("U")
            elif(i == "T"):
                rna.append("A")
            elif(i == "C"):
                rna.append("G")
            elif(i == "G"):
                rna.append("C")
            else:
                pass
        return "".join(rna)   

class Replicacao(DNA):
    def __init__(self,sequencia):
        super().__init__(sequencia)  

    def polimerase(self):
        dna = []
        for i in self._sequencia:
            if(i == "A"):
                dna.append("T")
            elif(i == "T"):
                dna.append("A")
            elif(i == "C"):
                dna.append("G")
            elif(i == "G"):
                dna.append("C")
            else:
                pass
        return "".join(dna)

    def fita_final(self):
        r = self.polimerase()
        return r + "-" + self._sequencia


class Traducao(Transcricao):
    def __init__(self,sequencia):
        super().__init__(sequencia)

    def __codon_aminoacido(self,codon):
        codon_tabela = {"UUU":"F","UUC":"F","UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L","AUC":"I","AUU":"I","AUA":"I","AUG":"M",
                 "GUU":"V","GUC":"V","GUA":"V","GUG":"V","UAC":"Y","UAU":"Y","UGU":"C","UGC":"C","UGG":"W","UCU":"S","UCC":"S","UCA":"S","UCG":"S"}
        return codon_tabela[codon]

    def traducao(self):
        rna = self.transcriptase()
        rib = []
        proteina = []
        for i in range(1,len(rna)+1):
            rib.append(rna[i-1])
            if(i%3==0 and i>0):
                codon = "".join(rib)
                if(codon == "UAA" or codon == "UAG" or codon == "UGA"):
                    break
                aminoacido = self.__codon_aminoacido(codon)
                proteina.append(aminoacido)
                rib = []
        return "-".join(proteina)


    
sequencia = "CACATGCATTATAAAATCCAA"
D = DNA(sequencia)
print(round(D.get_estabilidade,2))
t = Transcricao(sequencia)
print(t.transcriptase())
print(t.get_dna)
#----------------------------------------------

print("Testando a REPLICAÇÃO")
print("-----------------------")
r = Replicacao(sequencia)
print(r.polimerase())
print(r.fita_final())
print(r.get_dna)
print("Testando a síntese de proteínas")
print("-----------------------")
trad = Traducao(sequencia)
print('Sequência de aminoácidos: {}'.format(trad.traducao()))
#print(DNA.codon_aminoacido())

