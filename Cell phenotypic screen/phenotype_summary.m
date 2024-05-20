a=unique(class_low,'stable')
b=cellfun(@(x) sum(ismember(class_low,x)),a1,'un',0)

a1=unique(class_high,'stable')
b1=cellfun(@(x) sum(ismember(class_high,x)),a1,'un',0)