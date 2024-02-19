@def title="Single transferable vote plus (STV+)"
@def tags=["electoral-systems"]
@def mintoclevel=1

This webpage is about an electoral system I have devised that I wanted to share. I call it single transferable vote plus (STV+), and as its name suggests it is an additional member system based on [single transferable vote](https://en.wikipedia.org/wiki/Single_transferable_vote).

\toc

# Single transferable vote
**Single transferable vote** (STV) is a voting system used in relatively few countries (namely Australia, Republic of Ireland, Northern Ireland, and Malta). It is used in elections in which each electorate is set to return multiple elected officials and voters are expected to assign preferences to candidates that are running in their seat.

## Mechanics
In it, a quota of votes a candidate needs to win a seat is calculated (usually either the [Droop quota](https://en.wikipedia.org/wiki/Droop_quota) or [Hare quota](https://en.wikipedia.org/wiki/Hare_quota)). If a candidate reaches this on first preference votes, they are immediately elected and any surplus of votes (i.e. votes that have above the required quota) they have is redistributed according to next preferences. If all the seats that need to be filled are not filled by this first step, the next step is to eliminate the least popular remaining candidate and redistribute their votes according to next preferences. If candidate(s) have reached quota, they are elected and any surplus of votes redistributed, and this process repeats until all the seats that need to be filled are filled.

## Elected officials per electorate
The exact number of elected officials per electorate differs, but is usually between 3 and 7 for lower houses and potentially more for upper houses (up to 22 for South Australian Legislative Council elections after a double dissolution). 

## Problems
The fewer the electorates, hence greater the number of elected officials per electorate, the more proportionally representative the system gets. The more elected officials per electorate, however, the larger the ballot paper and the more candidates that need to be numbered by voters in order for the results to be meaningful. 

One way around this issue is to get voters to vote for parties instead of candidates directly and let the parties decide who fills the seats their seat. This, however, eliminates one of the main arguments in favour of single transferable vote, namely that it is non-partisan as it allows voters to vote for candidates directly and allows independents to run. Single transferable vote plus is intended to overcome these problems.

# Single transferable vote plus
Single transferable vote plus comprises two components:

* The single transferable vote component, where voters vote for candidates running in their multi-member electorates.
* The plus component, where voters assign preferences to parties and the proportion of the vote they get after preferences is used to decide the proportion of the seats they get. To proportion out the chamber, additional officials (who fill "overhang seats") are elected and they are the most popular losers from their party in the electorate races. 

The popularity of losing candidates is decided by first eliminating all winners from the race, then eliminating all the candidates from other parties and transfering all votes according to next preferences to the remaining candidates. Then eliminate the least popular excess candidates from each electorate and transfer votes according to next preferences (let $m$ be the number of overhang seats to be filled for that party, then every electorate will have at most $m$ candidates left after this step). Then whichever candidates in the electorate races have the most votes fill the party's overhang seats. 