# PBAT e MPFI-24: o que cada um pede de uma validação

## Sumário

- 10 (linha 15)
- 94 (linha 99)

Índice: PBAT · MPFI-24

O PBAT é o instrumento do mestrado planejado; o MPFI-24 fica aqui para análises
avulsas. Os dois pedem desenhos diferentes, e a diferença começa no que os
próprios autores dizem que o instrumento é. Fontes conferidas em 10/10/2026;
páginas citadas são do PDF do artigo.

## PBAT

O que o artigo original diz (Ciarrochi, Sahdra, Hofmann & Hayes, 2022):

- 18 itens, dois por linha da tabela do EEMM, um de valência positiva e um
  negativa: variação (2), seleção (12, em afeto, cognição, atenção, conexão
  social, motivação/autonomia e comportamento/competência), saúde física (2) e
  retenção (2) (p. 203).
- Resposta numa escala de 0 a 100 (discordo totalmente a concordo
  totalmente), escolhida no lugar de 5 ou 7 pontos para captar variabilidade
  intraindividual (p. 203).
- Amostra online transversal dos EUA, N = 598 (p. 202).
- Os autores afirmam que o PBAT não é uma escala, e sim um banco de itens para
  uso idionômico (p. 202), a ser usado como "gaveta de ferramentas", com mais
  ou menos itens conforme o caso (p. 208).
- Não houve análise fatorial nem de confiabilidade: os autores as consideram
  pouco adequadas para avaliar itens isolados, e a seleção de itens usou
  Boruta (pp. 202-204). As correlações médias foram .38 entre itens positivos,
  .49 entre negativos e perto de -.05 entre valências, o que os autores leem
  como valências distintas, sem testar modelo fatorial formal (pp. 204-205).
- Limitação declarada: o estudo transversal não trata o erro ergódico (p. 209).

Estudo intraindividual (Sanford et al., 2022): 50 pessoas responderam PBAT e
desfechos clínicos por 35 dias em amostragem de experiência, com pelo menos 60
ocasiões por pessoa, analisados com GIMME. As respostas violaram os
pressupostos ergódicos, não houve padrão geral de grupo, e a maior parte das
relações entre processo e desfecho foi idiográfica (resumo do artigo).

Adaptações publicadas que a skill localizou (texto lido só o da romena):

- Suécia: Larsson e Sundström (2024), abordagem centrada em Boruta (pelo
  título; texto não lido).
- Romênia: Popa et al. (2025), N = 637 (485 não clínicos, 152 clínicos).
  Redes EBICglasso sobre Pearson, CS de força .67 na rede não clínica e .20 na
  clínica, NCT entre as duas e, só na amostra não clínica, Boruta contra
  PHQ-9 e GAD-7. Os autores consideram o alfa conceitualmente inadequado
  porque os itens captam processos distintos, não escalas homogêneas. Houve
  piloto da tradução com 10 estudantes de medicina. Limitações declaradas:
  tradução sem comitê de especialistas, sem teste-reteste, grupo clínico
  pequeno e heterogêneo. O texto cita ainda adaptações polonesa e espanhola,
  não conferidas aqui.
- Coreia: Jo et al. (2025), abordagem de aprendizado de máquina (pelo título;
  texto não lido).

O que isso implica para a validação brasileira (leitura desta skill, para
discutir na orientação, não decisão):

1. AFC do PBAT como escala de dois fatores contraria o enquadramento dos
   autores (Ciarrochi et al., 2022, pp. 202 e 208) e o desenho da adaptação
   romena (Popa et al., 2025); as demais adaptações não foram lidas. Se a
   banca pedir estrutura, ela entra como descrição da organização por valência
   entre pessoas (a hipótese 1 do artigo original), não como validação de
   escore somado.
2. Um desenho coerente com o instrumento junta: adaptação transcultural pelas
   diretrizes da International Test Commission (2018), com comitê de juízes (a
   lacuna admitida pelos romenos) e pré-teste cognitivo;
   estrutura entre pessoas descritiva (correlações por valência, rede com CS,
   EGA e bootEGA, em Pearson porque a escala é 0-100); relevância dos itens
   para desfechos clínicos com instrumentos validados no Brasil (Boruta e
   correlações); comparação entre grupos com NCT e discussão de poder; e, para
   a evidência mais forte, um componente intraindividual em amostragem de
   experiência, que é o uso para o qual o PBAT foi desenhado.
3. Formato 0-100: Pearson e MLR, não policórica.
4. Confiabilidade: alfa e ômega de escore total não se aplicam a banco de
   itens. Rota da v1.4: estabilidade teste-reteste por item com
   `confiabilidade_teste_reteste()` (ICC(A,1), EPM e mudança mínima detectável)
   e, no componente de amostragem de experiência, `confiabilidade_ema()`, que
   separa a confiabilidade entre pessoas (RkF) da confiabilidade da mudança dentro
   da pessoa (Rc; Shrout & Lane, 2012; Cranford et al., 2006). Rc é a que sustenta
   análise intraindividual; com ela baixa, rede temporal e mudança ficam atenuadas.
   Petersen (2024, cap. 5) trata a teoria da generalizabilidade para medidas
   repetidas, e Martinková e Hladká (2023, cap. 4) os componentes de variância no
   R, os dois na pasta Psicometria-R.
5. Redes intraindividuais: o PBAT pede rede temporal e contemporânea idiográfica
   (N = 1) e multinível. Desde a v1.7, `rede_temporal()` (`references/redes-ega.md`,
   seção 10c) faz as duas pelos capítulos 9 e 10 de Isvoranu et al. (2022): com uma
   pessoa, graphicalVAR (PDC temporal e PCC contemporânea, EBIC); com várias, mlVAR
   (efeitos fixos temporais com p, DP dos aleatórios, redes contemporânea e entre
   pessoas). Destendencia por padrão e avisa quando há menos de 100 pontos (N = 1) ou
   menos de 20 por pessoa (multinível). Ficam fora: redes com variação no tempo (`mgm`,
   cap. 11), GIMME e o GVAR de painel do psychonetrics; esses seguem como rota fora
   do escopo testado. Antes da rede, `confiabilidade_ema()`: com Rc baixa as arestas
   temporais ficam atenuadas.

## MPFI-24

- MPFI completo (Rolffs, Rogge & Wilson, 2018): 60 itens em 12 subescalas de 5
  itens, 6 de flexibilidade e 6 de inflexibilidade do hexaflex, com itens
  selecionados por TRI. A pontuação divulgada pela ACBS atribui de 1 a 6
  pontos às respostas; conferir âncoras e período de referência no artigo.
- MPFI-24 (Grégoire et al., 2020): versões em inglês e francês. Pelo resumo, a
  AFC sustentou dois fatores de segunda ordem (flexibilidade e inflexibilidade)
  sobre seis fatores de primeira ordem cada, com amostras da população geral
  dos EUA (N = 2.668), de universitários francófonos do Canadá e da França
  (N = 728) e de trabalhadores francófonos de cinco países (N = 450). O
  Crossref registra uma errata do artigo: ler antes de citar números.

O que isso implica:

1. Modelo de segunda ordem com 12 fatores de 2 itens. Os evals da skill rodam
   exatamente esse formato (24 itens, 6 categorias, n = 1000) e recuperam as
   cargas de segunda ordem e a correlação entre os fatores gerais.
2. Seis categorias: com n acima de 200 e sem assimetria forte, WLSMV principal
   e MLR como sensibilidade (com 6 ou 7 categorias os dois são aceitáveis em
   muitas condições; Rhemtulla et al., 2012); com amostra pequena e limiares
   simétricos, MLR principal; com assimetria forte, WLSMV.
3. Modelos concorrentes a comparar, decididos antes da análise: 12 fatores
   correlacionados, dois fatores de segunda ordem, e o que a teoria do
   hexaflex justificar além disso.
4. Fatores de 2 itens têm confiabilidade limitada e AVE e HTMT pouco
   informativos; o relato de confiabilidade se apoia nos compostos gerais,
   com `confiabilidade_compostos()` nos dois modelos (ômega hierárquico no de
   segunda ordem, ômega total no de fatores correlacionados; ver
   `references/afc-ordinal.md`). v1.4: ajustar também o bifator (dois fatores
   gerais pedem dois bifatores, um por polo, ou um bifator com fator geral de
   flexibilidade psicológica, decidido pela teoria antes) e ler ωH, ωHS, ECV e PUC
   com `indices_bifator()` (Reise et al., 2018, cap. 22 do Wiley Handbook). Com
   subescalas de 2 itens, ωHS baixo é o esperado e diz que o escore da subescala
   toma emprestada quase toda a confiabilidade do fator geral.
5. Invariância (sexo, clínico e não clínico, faixa etária) com a sequência de
   Wu e Estabrook (2016) no modelo de fatores correlacionados de primeira
   ordem (`invariancia()` recusa o modelo de segunda ordem; o nível de segunda
   ordem pede passos montados à mão, ver `references/invariancia.md`), e TRI
   por subescala se houver pergunta sobre itens.
6. Antes de propor, conferir com a pesquisa-academica se já existe adaptação
   brasileira do MPFI ou do MPFI-24. A skill não conferiu.
