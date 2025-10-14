import os, smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText


def send_email(cfg, subject, text_body, html_body):
    em = cfg["output"]["email"]
    msg = MIMEMultipart("alternative")
    msg["Subject"] = subject
    msg["From"] = em["from_addr"]
    msg["To"] = ", ".join(em["to_addrs"])
    msg.attach(MIMEText(text_body, "plain", "utf-8"))
    msg.attach(MIMEText(html_body, "html", "utf-8"))

    password = os.environ.get(em["password_env"], "")
    with smtplib.SMTP(em["smtp_host"], em["smtp_port"]) as s:
        if em.get("use_starttls", True):
            s.starttls()
        s.login(em["username"], password)
        s.sendmail(em["from_addr"], em["to_addrs"], msg.as_string())
