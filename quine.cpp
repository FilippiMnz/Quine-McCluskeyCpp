#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <chrono>

struct Termo {
    std::string bits;
    std::vector<int> mintermos;
    bool coberto = false;

    Termo() {}
    Termo(int mintermo, int num_vars) {
        mintermos.push_back(mintermo);
        bits = "";
        for (int i = num_vars - 1; i >= 0; --i) {
            bits += ((mintermo >> i) & 1) ? '1' : '0';
        }
    }

    Termo(std::string b, std::vector<int> m) : bits(b), mintermos(m) {
        std::sort(mintermos.begin(), mintermos.end());
        mintermos.erase(std::unique(mintermos.begin(), mintermos.end()), mintermos.end());
    }

    bool operator==(const Termo& outro) const {
        return bits == outro.bits && mintermos == outro.mintermos;
    }
};

int contarUns(const std::string& s) {
    return std::count(s.begin(), s.end(), '1');
}

Termo combinarTermos(const Termo& t1, const Termo& t2) {
    int diferencas = 0;
    int indice_diferenca = -1;
    for (size_t i = 0; i < t1.bits.length(); ++i) {
        if (t1.bits[i] != t2.bits[i]) {
            diferencas++;
            indice_diferenca = i;
        }
    }

    if (diferencas == 1) {
        std::string novos_bits = t1.bits;
        novos_bits[indice_diferenca] = '-';
        std::vector<int> mintermos_combinados = t1.mintermos;
        mintermos_combinados.insert(mintermos_combinados.end(), t2.mintermos.begin(), t2.mintermos.end());
        return Termo(novos_bits, mintermos_combinados);
    }
    return Termo();
}


std::string paraAlgebrico(const Termo& termo, int num_vars) {
    std::string resultado = "(";
    for (size_t i = 0; i < termo.bits.length(); ++i) {
        if (termo.bits[i] == '1') {
            resultado += "v" + std::to_string(i + 1);
        } else if (termo.bits[i] == '0') {
            resultado += "!" + std::string("v") + std::to_string(i + 1);
        }
        if (termo.bits[i] != '-' && i < termo.bits.length() - 1) {
            resultado += " "; 
        }
    }
    resultado += ")";
    return resultado == "()" ? "1" : resultado;
}

bool analisarPLA(const std::string& nome_ficheiro, int& num_vars, std::vector<int>& mintermos) {
    std::ifstream ficheiro(nome_ficheiro);
    if (!ficheiro.is_open()) {
        std::cerr << "Erro: Não foi possível abrir o ficheiro " << nome_ficheiro << std::endl;
        return false;
    }

    std::string linha;
    while (std::getline(ficheiro, linha)) {
        std::stringstream ss(linha);
        std::string palavra_chave;
        ss >> palavra_chave;

        if (palavra_chave == ".i") {
            ss >> num_vars;
        } else if (palavra_chave == ".o" || palavra_chave == ".p" || palavra_chave == ".ilb" || palavra_chave == ".ob") {
            continue;
        } else if (palavra_chave == ".e") {
            break;
        } else {
            std::string bits_termo = palavra_chave;
            std::string saida;
            ss >> saida;

            if (saida == "1") {
                try {
                    int mintermo = std::stoi(bits_termo, nullptr, 2);
                    mintermos.push_back(mintermo);
                } catch (...) {
                }
            }
        }
    }
    ficheiro.close();
    std::sort(mintermos.begin(), mintermos.end());
    mintermos.erase(std::unique(mintermos.begin(), mintermos.end()), mintermos.end());
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <ficheiro.pla>" << std::endl;
        return 1;
    }

    int num_variaveis = 0;
    std::vector<int> vetor_mintermos;

    if (!analisarPLA(argv[1], num_variaveis, vetor_mintermos) || num_variaveis == 0) {
        std::cerr << "Erro no processamento do ficheiro PLA." << std::endl;
        return 1;
    }

    std::set<int> conjunto_mintermos(vetor_mintermos.begin(), vetor_mintermos.end());

    std::cout << "Número de variáveis: " << num_variaveis << std::endl;
    std::cout << "Mintermos: ";
    for (int m : conjunto_mintermos) std::cout << m << " ";
    std::cout << "\n---\n";

    auto inicio = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<Termo>> grupos(num_variaveis + 1);
    for (int m : conjunto_mintermos) {
        Termo t(m, num_variaveis);
        grupos[contarUns(t.bits)].push_back(t);
    }

    std::vector<Termo> implicantes_primos;
    std::vector<std::vector<Termo>> grupos_atuais = grupos;

    while (true) {
        for (auto& grupo : grupos_atuais) {
            for (auto& termo : grupo) termo.coberto = false; 
        }

        std::vector<std::vector<Termo>> proximos_grupos(num_variaveis + 1);
        bool mudou = false;

        for (size_t i = 0; i < grupos_atuais.size() - 1; ++i) {
            for (Termo& t1 : grupos_atuais[i]) {
                for (Termo& t2 : grupos_atuais[i + 1]) {
                    Termo combinado = combinarTermos(t1, t2);
                    if (!combinado.bits.empty()) {
                        t1.coberto = true;
                        t2.coberto = true;
                        proximos_grupos[contarUns(combinado.bits)].push_back(combinado);
                        mudou = true;
                    }
                }
            }
        }

        for (const auto& grupo : grupos_atuais) {
            for (const Termo& termo : grupo) {
                if (!termo.coberto) {
                    implicantes_primos.push_back(termo);
                }
            }
        }

        if (!mudou) break;

        for (auto& grupo : proximos_grupos) {
            std::sort(grupo.begin(), grupo.end(), [](const Termo& a, const Termo& b){ return a.bits < b.bits; });
            grupo.erase(std::unique(grupo.begin(), grupo.end()), grupo.end());
        }
        grupos_atuais = proximos_grupos;
    }

    std::sort(implicantes_primos.begin(), implicantes_primos.end(), [](const Termo& a, const Termo& b){ return a.bits < b.bits; });
    implicantes_primos.erase(std::unique(implicantes_primos.begin(), implicantes_primos.end()), implicantes_primos.end());

    std::cout << "Implicantes Primos:\n";
    for (const auto& ip : implicantes_primos) {
        std::cout << "  " << paraAlgebrico(ip, num_variaveis) << " (" << ip.bits << ") cobre m";
        for (int m : ip.mintermos) std::cout << m << ",";
        std::cout << "\b \n";
    }
    std::cout << "---\n";

    std::map<int, std::vector<int>> tabela;
    for (size_t i = 0; i < implicantes_primos.size(); ++i) {
        for (int m : implicantes_primos[i].mintermos) {
            tabela[m].push_back(i);
        }
    }

    std::vector<Termo> solucao;
    std::set<int> mintermos_cobertos;
    std::vector<bool> ip_utilizado(implicantes_primos.size(), false);

    for (const auto& par : tabela) {
        if (par.second.size() == 1) {
            int indice_ip = par.second[0];
            if (!ip_utilizado[indice_ip]) {
                solucao.push_back(implicantes_primos[indice_ip]);
                ip_utilizado[indice_ip] = true;
                for (int m : implicantes_primos[indice_ip].mintermos) {
                    mintermos_cobertos.insert(m);
                }
            }
        }
    }

    std::cout << "Implicantes Essenciais:\n";
    for (const auto& termo : solucao) {
        std::cout << "  " << paraAlgebrico(termo, num_variaveis) << "\n";
    }
    std::cout << "---\n";

    std::set<int> mintermos_restantes;
    for (int m : conjunto_mintermos) {
        if (!mintermos_cobertos.count(m)) {
            mintermos_restantes.insert(m);
        }
    }

    while (!mintermos_restantes.empty()) {
        int melhor_indice_ip = -1;
        int max_cobertos = 0;

        for (size_t i = 0; i < implicantes_primos.size(); ++i) {
            if (ip_utilizado[i]) continue;

            int cobertos_atualmente = 0;
            for (int m : implicantes_primos[i].mintermos) {
                if (mintermos_restantes.count(m)) {
                    cobertos_atualmente++;
                }
            }

            if (cobertos_atualmente > max_cobertos) {
                max_cobertos = cobertos_atualmente;
                melhor_indice_ip = i;
            }
        }

        if (melhor_indice_ip != -1) {
            solucao.push_back(implicantes_primos[melhor_indice_ip]);
            ip_utilizado[melhor_indice_ip] = true;
            for (int m : implicantes_primos[melhor_indice_ip].mintermos) {
                mintermos_restantes.erase(m);
            }
        } else {
            std::cout << "Não foi possível cobrir todos os mintermos.\n";
            break;
        }
    }

    std::cout << "Função Minimizada Final:\nF = ";
    for (size_t i = 0; i < solucao.size(); ++i) {
        std::cout << paraAlgebrico(solucao[i], num_variaveis);
        if (i < solucao.size() - 1) std::cout << " + ";
    }
    std::cout << "\n";

    auto fim = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duracao = fim - inicio;
    std::cout << "Tempo de minimização: " << duracao.count() << " ms\n";

    return 0;
}

