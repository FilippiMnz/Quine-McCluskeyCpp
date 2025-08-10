import sys
from collections import defaultdict

class Termo:
    def __init__(self, bits="", mintermos=None):
        self.bits = bits
        self.mintermos = sorted(mintermos) if mintermos else []
        self.coberto = False

    def __eq__(self, other):
        return self.bits == other.bits and self.mintermos == other.mintermos

    def __lt__(self, other):
        return self.bits < other.bits

def contar_uns(bits):
    return bits.count('1')

def combinar_termos(t1, t2):
    diferencas = 0
    indice_diferenca = -1
    for i in range(len(t1.bits)):
        if t1.bits[i] != t2.bits[i]:
            diferencas += 1
            indice_diferenca = i
    if diferencas == 1:
        novos_bits = list(t1.bits)
        novos_bits[indice_diferenca] = '-'
        mintermos_combinados = sorted(set(t1.mintermos + t2.mintermos))
        return Termo("".join(novos_bits), mintermos_combinados)
    return None

def para_algebrico(termo):
    resultado = []
    for i, bit in enumerate(termo.bits):
        if bit == '1':
            resultado.append(f"v{i+1}")
        elif bit == '0':
            resultado.append(f"!v{i+1}")
    if not resultado:
        return "1"
    return "(" + " ".join(resultado) + ")"

def analisar_pla(nome_arquivo):
    num_vars = 0
    mintermos = []
    try:
        with open(nome_arquivo, 'r') as f:
            for linha in f:
                linha = linha.strip()
                if not linha or linha.startswith('#'):
                    continue
                partes = linha.split()
                if partes[0] == '.i':
                    num_vars = int(partes[1])
                elif partes[0] in ['.o', '.p', '.ilb', '.ob']:
                    continue
                elif partes[0] == '.e':
                    break
                else:
                    if len(partes) >= 2 and partes[1] == '1':
                        bits = partes[0]
                        try:
                            mintermo = int(bits, 2)
                            mintermos.append(mintermo)
                        except:
                            pass
        mintermos = sorted(set(mintermos))
        return num_vars, mintermos
    except FileNotFoundError:
        print(f"Erro: arquivo '{nome_arquivo}' não encontrado.")
        return 0, []

def quine_mccluskey(num_vars, mintermos):
    grupos = [[] for _ in range(num_vars+1)]
    for m in mintermos:
        bits = format(m, f"0{num_vars}b")
        grupos[contar_uns(bits)].append(Termo(bits, [m]))

    implicantes_primos = []
    grupos_atuais = grupos

    while True:
        for grupo in grupos_atuais:
            for termo in grupo:
                termo.coberto = False

        proximos_grupos = [[] for _ in range(num_vars+1)]
        mudou = False

        for i in range(len(grupos_atuais)-1):
            for t1 in grupos_atuais[i]:
                for t2 in grupos_atuais[i+1]:
                    combinado = combinar_termos(t1, t2)
                    if combinado:
                        t1.coberto = True
                        t2.coberto = True
                        if combinado not in proximos_grupos[contar_uns(combinado.bits)]:
                            proximos_grupos[contar_uns(combinado.bits)].append(combinado)
                        mudou = True

        for grupo in grupos_atuais:
            for termo in grupo:
                if not termo.coberto and termo not in implicantes_primos:
                    implicantes_primos.append(termo)

        if not mudou:
            break

        grupos_atuais = proximos_grupos

    implicantes_primos.sort()

    tabela = defaultdict(list)
    for i, ip in enumerate(implicantes_primos):
        for m in ip.mintermos:
            tabela[m].append(i)

    solucao = []
    mintermos_cobertos = set()
    ip_utilizado = [False] * len(implicantes_primos)
    for m, ips in tabela.items():
        if len(ips) == 1:
            idx = ips[0]
            if not ip_utilizado[idx]:
                solucao.append(implicantes_primos[idx])
                ip_utilizado[idx] = True
                mintermos_cobertos.update(implicantes_primos[idx].mintermos)

    mintermos_restantes = set(mintermos) - mintermos_cobertos
    while mintermos_restantes:
        melhor_idx = -1
        max_cobertos = 0
        for i, ip in enumerate(implicantes_primos):
            if ip_utilizado[i]:
                continue
            cobertos = len(mintermos_restantes.intersection(ip.mintermos))
            if cobertos > max_cobertos:
                max_cobertos = cobertos
                melhor_idx = i
        if melhor_idx == -1:
            print("Não foi possível cobrir todos os mintermos.")
            break
        solucao.append(implicantes_primos[melhor_idx])
        ip_utilizado[melhor_idx] = True
        mintermos_restantes -= set(implicantes_primos[melhor_idx].mintermos)

    return [para_algebrico(t) for t in solucao]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Uso: python {sys.argv[0]} <arquivo.pla>")
        sys.exit(1)

    arquivo = sys.argv[1]
    num_vars, mintermos = analisar_pla(arquivo)
    if num_vars == 0 or not mintermos:
        print("Erro ao processar o arquivo PLA.")
        sys.exit(1)

    resultado = quine_mccluskey(num_vars, mintermos)
    print("Função Minimizada Final:")
    print("F = " + " + ".join(resultado))

