# coding: utf-8

"""
Small module to handle email notifications
Currently only basic functionality
"""

import re as re
import smtplib as smtplib
from email.mime.text import MIMEText


# make it simple here; we all know that specifying a regexp
# that matches all valid email addresses is not possible ;-)
__RE_EMAIL_ADDRESS__ = "^[\w\.-]+@[\w\.-]+\.\w{2,4}$"


def send_notification(address, subname, exitcode, start, end, log):
    """
    :param address: email address
    :param subname: name of executed subcommand
    :param exitcode: exit code of run
    :param start: date/time
    :param end: date/time
    :param log: buffer holding logging info
     :type: io.StringIO
    :return:
    """
    content = 'Hello {user},\nyour CREEPIEST run {name} finished.\n' \
              'Start: {starttime}\nEnd: {endtime}\nLog:\n{logoutput}\n\n' \
              'Kind regards,\nthe CREEPIEST framework'
    kwargs = {'user': address.split('@')[0], 'starttime': start, 'endtime': end, 'name': subname}
    logoutput = '<not requested>'
    if log:
        logoutput = log.getvalue()
    kwargs['logoutput'] = logoutput
    mobj = re.search(__RE_EMAIL_ADDRESS__, address)
    if mobj is not None:
        try:
            msg = MIMEText(content.format(**kwargs))
            msg['Subject'] = 'CREEPIEST run exit code: {}'.format(exitcode)
            msg['From'] = 'DoNotReply@creepiest-tool.org'
            msg['To'] = address
            s = smtplib.SMTP('localhost')
            s.send_message(msg)
            s.quit()
        except Exception as e:
            # whatever happens, do not risk disgraceful shutdown at this point
            # should it write to stderr?
            pass
    return
