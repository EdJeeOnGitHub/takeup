# Codebook: endline item-meaning free-text responses

Version 1.0 — locked before any coding pass. Any revision requires re-coding every
string in both passes.

## Task

Each row of `responses-dedup.csv` is one distinct free-text answer to the endline
question asking what it means when someone has the item (`item` column: bracelet, ink,
or calendar). Assign one or more of the five category labels below to every
`response_id`. Responses are written by field enumerators transcribing spoken answers;
expect Kenyan-English phrasing, misspellings ("derwomed", "medicen", "trested"), and
third- or first-person framing. Code the *content*, not the spelling or person.

Output format: CSV with columns `response_id,M,P,C,V,O`, one row per input row, each
label column 0 or 1, at least one label per row.

## Categories

### M — motivation / character attribution
The response attributes a trait, motive, or disposition to the person with the item:
they care about their health, are responsible, health-conscious, hygienic, set an
example, care about the community, are serious about prevention.

- "One is taking care of his Health" -> M (also P if participation is stated)
- "It shows they care about their health" -> M + P
- "A responsible person" -> M
- "Sets an example to others" -> M

The key test: does the answer tell you something about *what kind of person* they are,
beyond the bare fact of having participated? "Protecting themselves from worms" states
a purpose of the action rather than a trait of the person: P only, not M, unless
wording marks it as characteristic care/conscientiousness ("cares about", "keen on",
"serious about", "takes care of").

### P — participation marker
The response says the item shows the person took deworming treatment / was treated /
went for deworming. The overwhelmingly common answer.

- "They have been dewormed" -> P
- "He was given the deworming tablet" -> P
- "They got derwomed" -> P
- "Means he went for treatment" -> P
- "Protects them from worms" / "prevention of worms" -> P (purpose of the action)

### C — campaign / program / verification function
The response describes an administrative or informational function: identifying who was
treated for the program, proof of treatment, avoiding double-dosing, association with
the campaign itself.

- "You cannot get drug twice" -> C
- "Proof that one is treated" -> C + P
- "It was being marked on everyone who received the drugs" -> C + P
- "From the deworming campaign" -> C

### V — private / decorative / functional item value
The response is about the item's own value to the holder: gift, appreciation, beauty,
decoration; for the calendar, telling dates or reminding about the next deworming.

- "Was a gift for taking the drug" -> V (+ P, since taking the drug is stated)
- "Aesthetic purposes" -> V
- "Its to remind on the next deworming date" -> V
- "It assists as to know dates and months of the year." -> V
- "Appreciation" -> V

### O — other / don't know / uncodable
Nothing codable above: "Don't know", "Doesn't know", bare numerals or item names
("10", "Bracelets"), off-topic or incomprehensible fragments, refusals.

- "1" -> O
- "Don't know" -> O
- "Was  not around does not" -> O

## Decision rules

1. Multi-label: assign every category whose content is present. "A present; & a sign
   that dewormed" -> V + P.
2. M requires an attribution about the person. Mentions of health as the *purpose* of
   deworming ("to be healthy", "kills worms") are P, not M. Wording such as "cares
   about", "takes care of", "concerned with", "keen on", "responsible", "example to
   others", "loves his health" is M.
3. P covers any statement that treatment happened, in any tense or person, however
   misspelled.
4. C is for identification/proof/no-repeat-dosing/campaign-branding functions; the
   plain statement "shows they were dewormed" is P alone (a marker), C only when the
   administrative or verification purpose is explicit.
5. O is exclusive: if O = 1, all other labels are 0.
6. Code only from the text given. Do not guess from the item type; a bare "Dewormed"
   under calendar is P exactly as under bracelet.
7. When genuinely torn between including M or not, do not assign M. (Conservative
   default keeps the primary outcome strict.)
