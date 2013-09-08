from lib.Monomer import Monomer

class TypeCys(Monomer):
    name = "Cys"
    interaction_table = { 'Cys': -5.4400,'Met': -4.9900,'Phe': -5.8000,'Ile': -5.5000,'Leu': -5.8300,'Val': -4.9600,'Trp': -4.9500,'Tyr': -4.1600,'Ala': -3.5700,'Gly': -3.1600,'Thr': -3.1100,'Ser': -2.8600,'Asn': -2.5900,'Gln': -2.8500,'Asp': -2.4100,'Glu': -2.2700,'His': -3.6000,'Arg': -2.5700,'Lys': -1.9500,'Pro': -3.0700 }

class TypeIle(Monomer):
    name = "Ile"
    interaction_table = { 
'Cys': -5.5000,'Met': -6.0200,'Phe': -6.8400,'Ile': -6.5400,'Leu': -7.0400,'Val': -6.0500,'Trp': -5.7800,'Tyr': -5.2500,'Ala': -4.5800,'Gly': -3.7800,'Thr': -4.0300,'Ser': -3.5200,'Asn': -3.2400,'Gln': -3.6700,'Asp': -3.1700,'Glu': -3.2700,'His': -4.1400,'Arg': -3.6300,'Lys': -3.0100,'Pro': -3.7600 }
    
class TypeSer(Monomer):
    name = "Ser"
    interaction_table = { 
'Cys': -2.8600,'Met': -3.0300,'Phe': -4.0200,'Ile': -3.5200,'Leu': -3.9200,'Val': -3.0500,'Trp': -2.9900,'Tyr': -2.7800,'Ala': -2.0100,'Gly': -1.8200,'Thr': -1.9600,'Ser': -1.6700,'Asn': -1.5800,'Gln': -1.4900,'Asp': -1.6300,'Glu': -1.4800,'His': -2.1100,'Arg': -1.6200,'Lys': -1.0500,'Pro': -1.5700 }
    
class TypeVal(Monomer):
    name = "Val"
    interaction_table = { 
'Cys': -4.9600,'Met': -5.3200,'Phe': -6.2900,'Ile': -6.0500,'Leu': -6.4800,'Val': -5.5200,'Trp': -5.1800,'Tyr': -4.6200,'Ala': -4.0400,'Gly': -3.3800,'Thr': -3.4600,'Ser': -3.0500,'Asn': -2.8300,'Gln': -3.0700,'Asp': -2.4800,'Glu': -2.6700,'His': -3.5800,'Arg': -3.0700,'Lys': -2.4900,'Pro': -3.3200 }
    
class TypeGly(Monomer):
    name = "Gly"
    interaction_table = { 
'Cys': -3.1600,'Met': -3.3900,'Phe': -4.1300,'Ile': -3.7800,'Leu': -4.1600,'Val': -3.3800,'Trp': -3.4200,'Tyr': -3.0100,'Ala': -2.3100,'Gly': -2.2400,'Thr': -2.0800,'Ser': -1.8200,'Asn': -1.7400,'Gln': -1.6600,'Asp': -1.5900,'Glu': -1.2200,'His': -2.1500,'Arg': -1.7200,'Lys': -1.1500,'Pro': -1.8700 }
    
class TypeGln(Monomer):
    name = "Gln"
    interaction_table = { 
'Cys': -2.8500,'Met': -3.3000,'Phe': -4.1000,'Ile': -3.6700,'Leu': -4.0400,'Val': -3.0700,'Trp': -3.1100,'Tyr': -2.9700,'Ala': -1.8900,'Gly': -1.6600,'Thr': -1.9000,'Ser': -1.4900,'Asn': -1.7100,'Gln': -1.5400,'Asp': -1.4600,'Glu': -1.4200,'His': -1.9800,'Arg': -1.8000,'Lys': -1.2900,'Pro': -1.7300 }
    
class TypePro(Monomer):
    name = "Pro"
    interaction_table = { 
'Cys': -3.0700,'Met': -3.4500,'Phe': -4.2500,'Ile': -3.7600,'Leu': -4.2000,'Val': -3.3200,'Trp': -3.7300,'Tyr': -3.1900,'Ala': -2.0300,'Gly': -1.8700,'Thr': -1.9000,'Ser': -1.5700,'Asn': -1.5300,'Gln': -1.7300,'Asp': -1.3300,'Glu': -1.2600,'His': -2.2500,'Arg': -1.7000,'Lys': -0.9700,'Pro': -1.7500 }
    
class TypeLys(Monomer):
    name = "Lys"
    interaction_table = { 
'Cys': -1.9500,'Met': -2.4800,'Phe': -3.3600,'Ile': -3.0100,'Leu': -3.3700,'Val': -2.4900,'Trp': -2.6900,'Tyr': -2.6000,'Ala': -1.3100,'Gly': -1.1500,'Thr': -1.3100,'Ser': -1.0500,'Asn': -1.2100,'Gln': -1.2900,'Asp': -1.6800,'Glu': -1.8000,'His': -1.3500,'Arg': -0.5900,'Lys': -0.1200,'Pro': -0.9700 }
    
class TypeThr(Monomer):
    name = "Thr"
    interaction_table = { 
'Cys': -3.1100,'Met': -3.5100,'Phe': -4.2800,'Ile': -4.0300,'Leu': -4.3400,'Val': -3.4600,'Trp': -3.2200,'Tyr': -3.0100,'Ala': -2.3200,'Gly': -2.0800,'Thr': -2.1200,'Ser': -1.9600,'Asn': -1.8800,'Gln': -1.9000,'Asp': -1.8000,'Glu': -1.7400,'His': -2.4200,'Arg': -1.9000,'Lys': -1.3100,'Pro': -1.9000 }

    
class TypePhe(Monomer):
    name = "Phe"
    interaction_table = { 
'Cys': -5.8000,'Met': -6.5600,'Phe': -7.2600,'Ile': -6.8400,'Leu': -7.2800,'Val': -6.2900,'Trp': -6.1600,'Tyr': -5.6600,'Ala': -4.8100,'Gly': -4.1300,'Thr': -4.2800,'Ser': -4.0200,'Asn': -3.7500,'Gln': -4.1000,'Asp': -3.4800,'Glu': -3.5600,'His': -4.7700,'Arg': -3.9800,'Lys': -3.3600,'Pro': -4.2500 }

class TypeAla(Monomer):
    name = "Ala"
    interaction_table = { 
'Cys': -3.5700,'Met': -3.9400,'Phe': -4.8100,'Ile': -4.5800,'Leu': -4.9100,'Val': -4.0400,'Trp': -3.8200,'Tyr': -3.3600,'Ala': -2.7200,'Gly': -2.3100,'Thr': -2.3200,'Ser': -2.0100,'Asn': -1.8400,'Gln': -1.8900,'Asp': -1.7000,'Glu': -1.5100,'His': -2.4100,'Arg': -1.8300,'Lys': -1.3100,'Pro': -2.0300 }

    
class TypeMet(Monomer):
    name = "Met"
    interaction_table = { 
'Cys': -4.9900,'Met': -5.4600,'Phe': -6.5600,'Ile': -6.0200,'Leu': -6.4100,'Val': -5.3200,'Trp': -5.5500,'Tyr': -4.9100,'Ala': -3.9400,'Gly': -3.3900,'Thr': -3.5100,'Ser': -3.0300,'Asn': -2.9500,'Gln': -3.3000,'Asp': -2.5700,'Glu': -2.8900,'His': -3.9800,'Arg': -3.1200,'Lys': -2.4800,'Pro': -3.4500 }

    
class TypeAsp(Monomer):
    name = "Asp"
    interaction_table = { 
'Cys': -2.4100,'Met': -2.5700,'Phe': -3.4800,'Ile': -3.1700,'Leu': -3.4000,'Val': -2.4800,'Trp': -2.8400,'Tyr': -2.7600,'Ala': -1.7000,'Gly': -1.5900,'Thr': -1.8000,'Ser': -1.6300,'Asn': -1.6800,'Gln': -1.4600,'Asp': -1.2100,'Glu': -1.0200,'His': -2.3200,'Arg': -2.2900,'Lys': -1.6800,'Pro': -1.3300 }

    
class TypeLeu(Monomer):
    name = "Leu"
    interaction_table = { 
'Cys': -5.8300,'Met': -6.4100,'Phe': -7.2800,'Ile': -7.0400,'Leu': -7.3700,'Val': -6.4800,'Trp': -6.1400,'Tyr': -5.6700,'Ala': -4.9100,'Gly': -4.1600,'Thr': -4.3400,'Ser': -3.9200,'Asn': -3.7400,'Gln': -4.0400,'Asp': -3.4000,'Glu': -3.5900,'His': -4.5400,'Arg': -4.0300,'Lys': -3.3700,'Pro': -4.2000 }

    
class TypeHis(Monomer):
    name = "His"
    interaction_table = { 
'Cys': -3.6000,'Met': -3.9800,'Phe': -4.7700,'Ile': -4.1400,'Leu': -4.5400,'Val': -3.5800,'Trp': -3.9800,'Tyr': -3.5200,'Ala': -2.4100,'Gly': -2.1500,'Thr': -2.4200,'Ser': -2.1100,'Asn': -2.0800,'Gln': -1.9800,'Asp': -2.3200,'Glu': -2.1500,'His': -3.0500,'Arg': -2.1600,'Lys': -1.3500,'Pro': -2.2500 }

    
class TypeArg(Monomer):
    name = "Arg"
    interaction_table = { 
'Cys': -2.5700,'Met': -3.1200,'Phe': -3.9800,'Ile': -3.6300,'Leu': -4.0300,'Val': -3.0700,'Trp': -3.4100,'Tyr': -3.1600,'Ala': -1.8300,'Gly': -1.7200,'Thr': -1.9000,'Ser': -1.6200,'Asn': -1.6400,'Gln': -1.8000,'Asp': -2.2900,'Glu': -2.2700,'His': -2.1600,'Arg': -1.5500,'Lys': -0.5900,'Pro': -1.7000 }

    
class TypeTrp(Monomer):
    name = "Trp"
    interaction_table = { 
'Cys': -4.9500,'Met': -5.5500,'Phe': -6.1600,'Ile': -5.7800,'Leu': -6.1400,'Val': -5.1800,'Trp': -5.0600,'Tyr': -4.6600,'Ala': -3.8200,'Gly': -3.4200,'Thr': -3.2200,'Ser': -2.9900,'Asn': -3.0700,'Gln': -3.1100,'Asp': -2.8400,'Glu': -2.9900,'His': -3.9800,'Arg': -3.4100,'Lys': -2.6900,'Pro': -3.7300 }

    
class TypeGlu(Monomer):
    name = "Glu"
    interaction_table = { 
'Cys': -2.2700,'Met': -2.8900,'Phe': -3.5600,'Ile': -3.2700,'Leu': -3.5900,'Val': -2.6700,'Trp': -2.9900,'Tyr': -2.7900,'Ala': -1.5100,'Gly': -1.2200,'Thr': -1.7400,'Ser': -1.4800,'Asn': -1.5100,'Gln': -1.4200,'Asp': -1.0200,'Glu': -0.9100,'His': -2.1500,'Arg': -2.2700,'Lys': -1.8000,'Pro': -1.2600 }

    
class TypeAsn(Monomer):
    name = "Asn"
    interaction_table = { 
'Cys': -2.5900,'Met': -2.9500,'Phe': -3.7500,'Ile': -3.2400,'Leu': -3.7400,'Val': -2.8300,'Trp': -3.0700,'Tyr': -2.7600,'Ala': -1.8400,'Gly': -1.7400,'Thr': -1.8800,'Ser': -1.5800,'Asn': -1.6800,'Gln': -1.7100,'Asp': -1.6800,'Glu': -1.5100,'His': -2.0800,'Arg': -1.6400,'Lys': -1.2100,'Pro': -1.5300 }

    
class TypeTyr(Monomer):
    name = "Tyr"
    interaction_table = { 
'Cys': -4.1600,'Met': -4.9100,'Phe': -5.6600,'Ile': -5.2500,'Leu': -5.6700,'Val': -4.6200,'Trp': -4.6600,'Tyr': -4.1700,'Ala': -3.3600,'Gly': -3.0100,'Thr': -3.0100,'Ser': -2.7800,'Asn': -2.7600,'Gln': -2.9700,'Asp': -2.7600,'Glu': -2.7900,'His': -3.5200,'Arg': -3.1600,'Lys': -2.6000,'Pro': -3.1900 }

