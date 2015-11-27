import os, sys
import requests

def count_pairs(sequence, secstruct):
    open = []
    pairs = {'gc': 0, 'gu': 0, 'ua': 0}
    for i in range(len(sequence)):
        if secstruct[i] == '(':
            open.append(sequence[i])
        elif secstruct[i] == ')':
            pair = open.pop() + sequence[i]
            if pair == "GC" or pair == "CG":
                pairs['gc'] += 1
            elif pair == "GU" or pair == "UG":
                pairs['gu'] += 1
            elif pair == "UA" or pair = "AU":
                pairs['ua'] += 1
    return pairs

def post_solutions(puzzleid, filename):
    filename = os.path.join(settings.PUZZLE_DIR, filename)
    for i, sol in enumerate(open(filename, 'r')):
        if sol.startswith('#'):
            continue
        post_solution(puzzleid, "eternabot solution %s" % i, sol.strip())
    return

def post_solution(puzzleid, title, sequence):
    fold = fold_utils.nupack_fold(sequence)
    pairs = count_pairs(sequence, fold[0][0])
    header = {'Content-Type': 'application/x-www-form-urlencoded'}
    login = {'type': 'login',
             'name': 'theeternabot',
             'pass': 'iamarobot',
             'workbranch': 'main'}
    solution = {'type': 'post_solution',
                'puznid': puzzleid,
                'title': title,
                'body': 'eternabot solution',
                'sequence': sequence,
                'energy': fold[1],
                'gc': pairs['gc'],
                'gu': pairs['gu'],
                'ua': pairs['ua'],
                'melt': 97,
                'pointsrank': 'false',
                'recommend-puzzle': 'true'}

    url = "http://www.eternagame.org"
    loginurl = "%s/login/" % url
    posturl = "%s/post/" % url
    with requests.Session() as s:
        r = s.post(loginurl, data=login, headers=header)
        r = s.post(posturl, data=solution, headers=header)
    return

def main():
    post_solutions(sys.argv[2], sys.argv[1])

if __name__ == "__main__":
    main()
